#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# Trim database to only include studies from target period and area.
# Correct species naming issues in data

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
  .dbPF <- file.path(.wd,'processed_data/mosey_mod_2023.db')
  
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

evt <- tbl(db,'event')  # %>%  collect()

#---- Perform analysis ----#

#-- Make a filtered table by study period

# Collect inidvidual table
indtb <- tbl(db, "individual") %>%  collect()

# Collect study table
stdtb <- tbl(db, "study") %>% collect()

# extract only relevant time periods
mod <- evt %>%
  filter((timestamp > !!periods$date[periods$cutpoint == "start_pre-ld_2019"] & 
            timestamp < !!periods$date[periods$cutpoint == "stop_2019"])
         | (timestamp > !!periods$date[periods$cutpoint == "start_pre-ld_2020"] & 
              timestamp < !!periods$date[periods$cutpoint == "stop_2020"])) %>%
  mutate(yr = strftime('%Y', timestamp),
         trt = case_when(
           timestamp >= !!periods$date[periods$cutpoint == "start_pre-ld_2019"] &
             timestamp < !!periods$date[periods$cutpoint == "start_ld_2019"] ~ "pre-ld_2019",
           timestamp >= !!periods$date[periods$cutpoint == "start_ld_2019"] &
             timestamp < !!periods$date[periods$cutpoint == "start_post-ld_2019"] ~ "ld_2019",
           timestamp >= !!periods$date[periods$cutpoint == "start_post-ld_2019"] &
             timestamp < !!periods$date[periods$cutpoint == "stop_2019"] ~ "post-ld_2019",
           timestamp >= !!periods$date[periods$cutpoint == "start_pre-ld_2020"] &
             timestamp < !!periods$date[periods$cutpoint == "start_ld_2020"] ~ "pre-ld_2020",
           timestamp >= !!periods$date[periods$cutpoint == "start_ld_2020"] &
             timestamp < !!periods$date[periods$cutpoint == "start_post-ld_2020"] ~ "ld_2020",
           timestamp >= !!periods$date[periods$cutpoint == "start_post-ld_2020"] &
             timestamp < !!periods$date[periods$cutpoint == "stop_2020"] ~ "post-ld_2020")
  ) %>% 
  collect()

mod2 <- mod %>% 
  left_join(indtb, by = "individual_id") %>% 
  # remove particular studies
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>% 
  filter(study_id != 1891172051) %>%
  filter(study_id != 1891403240) %>% 
  mutate(ind_f = as.factor(individual_id),
         species = taxon_canonical_name)%>%  # create factor version of ind for REs)
  mutate(species = case_when( # correct species names
    study_id == 1442516400 ~ "Anser caerulescens",
    study_id == 1233029719 ~ "Odocoileus virginianus",
    study_id == 1631574074 ~ "Ursus americanus",
    study_id == 1418296656 ~ "Numenius americanus",
    study_id == 474651680  ~ "Odocoileus virginianus",
    study_id == 1044238185 ~ "Alces alces",
    study_id == 2548691779 ~ "Odocoileus hemionus", # !!! CORRECT IN WORKFLOW !!! #
    study_id == 2575515057 ~ "Cervus elaphus",  # !!! CORRECT IN WORKFLOW !!! #
    TRUE ~ species
  ))%>% 
  mutate(species = case_when(
    species == "Chen caerulescens" ~ "Anser caerulescens",
    TRUE ~ species
  ))

# write results back to db
dbWriteTable(conn = db, name = "event_trim", value = mod2, append = F, overwrite = T)

# Correct the individual table
ind_out <- indtb %>% 
  semi_join(mod2, by = "individual_id")   # only retain inds contained in cleaned event table

# write table back to db
dbWriteTable(conn = db, name = "individual_trim", value = ind_out, append = FALSE, overwrite = T)


# Correct the study table 
std_out <- stdtb %>% 
  semi_join(ind_out, by = "study_id") # only retain studies contained in cleaned individual table

# write table back to db
dbWriteTable(conn = db, name = "study_trim", value = std_out, append = FALSE, overwrite = T)


#---- Finalize script ----#

# disconnect from db
dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))