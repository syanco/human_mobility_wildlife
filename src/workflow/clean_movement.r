#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# This script cleans data for the COVID-19 Animal Movement Project


'
Prep and clean data for the COVID-19 Animal Movement Project
See project documentation for details about anticipated directory structure.
Expects db to be of format mosey_db, with a table named "event"

Usage:
clean_movement.r [--db=<db>] 
clean_movement.r (-h | --help)

Control files:
  analysis/ctfs/dates.csv
  
Conda Environment: spatial

Options:
-h --help     Show this screen.
-v --version     Show version.
-d --db=<db> Path to movement database. Defaults to <wd>/data/move.db
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '/Users/scottyanco/Documents/covid-19_movement/'
  .dbPF <- file.path(.wd,'processed_data/mosey_mod_2023_stash.db')
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$db)
  
}

#---- Initialize Environment ----#

# Start time
t0 <- Sys.time()

# Run startup
source(file.path(.wd,'analysis/src/startup.r'))

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(sf)
    library(geosphere)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

`%notin%` <- Negate(`%in%`)

#---- Load control files ----#
periods <- read_csv(file.path(.wd,'analysis/ctfs/dates.csv'),
                    col_types=list("date" = col_date(format = "%m/%d/%Y"))) 

#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))



# dbBegin(db)

#---- Perform analysis ----#

#-- Clean Outliers

evt_trm <- tbl(db,'event_final') %>%  collect()
ind_trm <- tbl(db, "individual_final") %>%  collect()
std_trm <- tbl(db, "study_final") %>%  collect()

# beepr::beep()

cnt <- evt_trm %>% 
  group_by(individual_id) %>% 
  summarise(n = n())

# get step lengths and turn angles across dataset
evt_calc <- evt_trm %>% 
  # left_join(cnt) %>% 
  # filter(n < 3) %>% 
  as.data.frame() %>% 
  # group run calcs per individual
  group_by(individual_id) %>% 
  arrange(timestamp) %>% 
  mutate(rn = row_number(),
         # lead_d = geometry[row_number()+1],
         # timestamp = parse_date_time(timestamp, "ymd_HMS"),
         lag_lon = dplyr::lag(lon, 1), #get pvs lon
         lag_lat = dplyr::lag(lat, 1), # get pvs lat
         sl = distGeo(cbind(lon,lat), cbind(lag_lon, lag_lat)), # step length (m)
         dt = as.numeric(difftime(timestamp, dplyr::lag(timestamp, 1)), units='secs'), # time diff (secs)
         v = sl/dt, # velocity (m/s)
         bearing = bearing(cbind(lon,lat), cbind(lag_lon, lag_lat)), # bearing since last point
         ta = 180-abs(180 - abs(bearing - dplyr::lag(bearing, 1)) %% 360), # turn angle from current bearing and last bearing
         wk = lubridate::week(timestamp))

# calculate qunatile-based cutoffs
cuts <- evt_calc %>% 
  # as.data.frame() %>% 
  group_by(individual_id) %>% 
  summarize(
    qta = quantile(ta, probs = 0.95, na.rm = T, names = F),
    qsl = quantile(sl, probs = 0.95, na.rm = T, names = F),
    qv  = quantile(v,  probs = 0.95, na.rm = T, names = F)
  ) 

# filter outliers
evt_out <- evt_calc %>% 
  # just join the cutpoints back to the dataset
  left_join(cuts) %>% 
  # conservative outlier thresh, must be past 95% quant for either sl and TA
  filter(
    v < qv & ta < qta,
    v < 25) %>% # remove observations faster than 25 m/s
  ungroup() 

# write table back to db
dbWriteTable(conn = db, name = "event_clean", value = evt_out, append = FALSE, overwrite = T)


#-- Update Individual and Event Tables

# Individual Table
ind_out <- ind_trm %>% 
  semi_join(evt_out, by = "individual_id")   # only retain inds contained in cleaned event table

# write table back to db
dbWriteTable(conn = db, name = "individual_clean", value = ind_out, append = FALSE, overwrite = T)

# Study Table
std_out <- std_trm %>% 
  semi_join(ind_out, by = "study_id") # only retain studies contained in cleaned individual table

# write table back to db
dbWriteTable(conn = db, name = "study_clean", value = std_out, append = FALSE, overwrite = T)



#-- Generate Sample Size Sumaries

# No of species
(no_sp <- evt_out %>% pull(species) %>% unique() %>% length())

# Inds per species
(inds_p_sp <- evt_out %>% group_by(species) %>% summarize(n = n_distinct(ind_f)))

# Fixes per species
(fix_p_sp <- evt_out %>% group_by(species) %>% summarize(n = n()))

# No of individuals
(no_inds <- evt_out %>% pull(ind_f) %>% unique() %>% length())

#-- Write our sample size summaries

message("Writing out sample size report...")

sink("out/sample_report_after_data_clean.txt")
print(glue("Description of sample sizes after all data cleaning, prior to estimating dBBMMs and MVNH Niche Breadths"))
print(glue("\n ##########################################################"))

print(glue("\n \n Total number of species in data: {no_sp}"))

print(glue("\n \n Total number of individuals in data: {no_inds}"))

print(glue("\n \n Individuals per species: \n {inds_p_sp}"))

print(glue("\n ##########################################################"))

print(glue("\n \n Total number of fixes in data: {sum(fix_p_sp$n)}"))

print(glue("\n \n Fixes per species: \n {fix_p_sp}"))

sink()

#---- Finalize script ----#

# disconnect from db
dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
