#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# This script checks for data completeness and trims data based on minimum data requirements


'
Filters data to minimums for sample sizes per week and number of individuals 
per species, reports out sample sizes

Usage:
filter_data_mins.r [--db=<db>] <wkmin> <minsp>
filter_data_mins.r (-h | --help)

Parameters
  wkmin:    Minimum sample size per week
  minsp:    Minimum individuals per species

Options:
-h --help     Show this screen.
-v --version     Show version.
-d --db=<db>  Path to movement database. Defaults to <wd>/data/move.db

' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  # .wd <- '~/project/covid-19_movement'
  .wd <- getwd()
  
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
  .wkmin <- 30
  .minsp <- 5
  
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$db)
  .wkmin <- as.numeric(ag$wkmin)
  .minsp <- as.numeric(ag$minsp)
}

#---- Initialize Environment ----#

# Start time
t0 <- Sys.time()

# Run startup
source(file.path(.wd,'src/startup.r'))

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(sf)
    library(geosphere)
    library(lubridate)
    library(tidyverse)
  }))

# Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)

# create operator that is opposite of %in%
`%notin%` <- Negate(`%in%`)

#---- Load control files ----#
periods <- read_csv(file.path(.wd,'ctfs/dates.csv'),
                    col_types=list("date" = col_date(format = "%m/%d/%Y"))) 

#---- Initialize database ----#
message("Connecting to database...")
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

#---- Perform analysis ----#

# dbListTables(db)

#-- Load Cleaned Data

message("Loading data...")
evt_cln <- tbl(db,'event_trim') %>%  collect()
# evt_fin <- tbl(db,'event_final') %>%  collect()
ind_cln <- tbl(db, "individual_trim") %>%  collect()
# ind_cln <- tbl(db, "individual_trim") %>%  collect()
std_cln <- tbl(db, "study_trim") %>%  collect()
# std_cln <- tbl(db, "study_trim") %>%  collect()

# beepr::beep()

# -- Remove incomplete cases

message(glue("Removing weeks without complete env annotations"))

evt_out2 <- evt_cln %>%
  mutate(timestamp2 = ymd_hms(timestamp),
         # timestamp_char <- as.character(timestamp2)
         wk = week(timestamp2)) %>% # add column for week component of timestamp
  drop_na(tmax, ndvi, elev) # drop any rows with NA values for annotation cols

# evt_comp <- evt_fix %>%
#   select(species, ind_f, tmax, ndvi, elev) %>%
#   filter(complete.cases(.))
# 
# evt_out <- evt_sp %>%
#   semi_join(evt_comp, by = "ind_f")

#-- Remove species with too few individuals

message(glue("Removing species below theshold of {.minsp} individuals..."))

# get ind count per species
# remove species that have less than the minimum # of individuals
sp_sum <- evt_out2 %>%
  group_by(species) %>%
  summarize(nind = length(unique(ind_f))) %>%
  filter(nind >= .minsp)

# retain events only for species that have the minimum # of individuals
evt_sp <- evt_out2 %>% 
  semi_join(sp_sum, by = "species")

#-- Remove weeks with too few fixes per week

message(glue("Removing weeks below threshold of {.wkmin} fixes per individual"))

# get fix count per week
fix_sum <- evt_sp %>% 
  group_by(ind_f, yr, wk) %>% 
  summarize(n= n()) %>% 
  filter((n >= .wkmin))

# retain events only for weeks that have the minimum number of fixes
evt_fix <- evt_sp %>% 
  semi_join(fix_sum, by = "ind_f")




#-- Write out individual table

message("Writing event table back to database...")

# write table back to db
dbWriteTable(conn = db, name = "event_final", value = evt_fix, append = FALSE, overwrite = T)

#-- Sync up study and individual tables

message("Updating and writing individual table back to database...")

# Individual Table
ind_out <- ind_cln %>% 
  semi_join(evt_fix, by = "individual_id")   # only retain inds contained in cleaned event table

# write table back to db
dbWriteTable(conn = db, name = "individual_final", value = ind_out, append = FALSE, overwrite = T)

message("Updating and writing study table back to database...")

# Study Table
std_out <- std_cln %>% 
  semi_join(ind_out, by = "study_id") # only retain studies contained in cleaned individual table

# write table back to db
dbWriteTable(conn = db, name = "study_final", value = std_out, append = FALSE, overwrite = T)


#-- Generate Sample Size Summaries

# No of species
(no_sp <- evt_fix %>% pull(species) %>% unique() %>% length())

# Inds per species
(inds_p_sp <- evt_fix %>% group_by(species) %>% 
    summarize(n = n_distinct(ind_f)) %>%
    janitor::adorn_totals("row") )

# Fixes per species
(fix_p_sp <- evt_cln %>% group_by(species) %>%
    summarize(n = n()) %>%
    janitor::adorn_totals("row"))

# No of individuals
(no_inds <- evt_fix %>% pull(ind_f) %>% unique() %>% length())

# No of studies
(no_stds <- evt_cln %>% 
    group_by(species) %>%
    summarize(n = n_distinct(study_id)) %>%
    janitor::adorn_totals("row") )

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

dbDisconnect(db)

message("Script Complete!")

