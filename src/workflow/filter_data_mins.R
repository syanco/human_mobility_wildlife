#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# This script checks for data completeness and trims data absed on minimum data requirements


'
Filters data to minimums, reports out sample sizes

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
  
  # .wd <- '~/projects/covid-19_movement'
  .wd <- '~/Documents/covid-19_movement/'
  
  .dbPF <- file.path(.wd,'processed_data/mosey_mod_2023.db')
  .wkmin <- 30
  .minsp <- 5
  
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$db)
  .wkmin <- ag$wkmin
  .minsp <- ag$minsp
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
message("Connecting to database...")
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

#---- Perform analysis ----#

# dbListTables(db)

#-- Load Cleaned Data

message("Loading data...")
evt_cln <- tbl(db,'event_clean') %>%  collect()
ind_cln <- tbl(db, "individual_clean") %>%  collect()
# ind_cln <- tbl(db, "individual_trim") %>%  collect()
std_cln <- tbl(db, "study_clean") %>%  collect()
# std_cln <- tbl(db, "study_trim") %>%  collect()

#-- Remove species with too few individuals

message(glue("Removing species below theshold of {.minsp} individuals..."))
# get ind count per species
sp_sum <- evt_cln %>%
  group_by(species) %>%
  summarize(nind = length(unique(ind_f))) %>%
  filter(nind >= .minsp) #require a minimum# of individuals

evt_sp <- evt_cln %>% 
  semi_join(sp_sum, by = "species")


#-- Remove weeks with too few fixes per week

message(glue("Removing weeks below threshold of {.wkmin} fixes per individual"))

fix_sum <- evt_sp %>% 
  group_by(ind_f, yr, wk) %>% 
  summarize(n= n()) %>% 
  filter((n >= .wkmin))

evt_out <- evt_sp %>% 
  semi_join(fix_sum, by = "ind_f")


#-- Write out individual table

message("Writing event table back to database...")

message("Writing out ")
# write table back to db
dbWriteTable(conn = db, name = "event_final", value = evt_out, append = FALSE, overwrite = T)

#-- Sync up study and individual tables

message("Updating and writing individual table back to database...")

# Individual Table
ind_out <- ind_cln %>% 
  semi_join(evt_out, by = "individual_id")   # only retain inds contained in cleaned event table

# write table back to db
dbWriteTable(conn = db, name = "individual_final", value = ind_out, append = FALSE, overwrite = T)

message("Updating and writing study table back to database...")

# Study Table
std_out <- std_cln %>% 
  semi_join(ind_out, by = "study_id") # only retain studies contained in cleaned individual table

# write table back to db
dbWriteTable(conn = db, name = "study_final", value = std_out, append = FALSE, overwrite = T)

#-- Write our sample size summaries

message("Writing out sample size report...")

sink("out/sample_report_after_data_clean.txt")
print(glue("Description of sample sizes after all data cleaning, prior to estimating dBBMMs and MVNH Niche Breadths"))
print(glue("\n ##########################################################"))

print(glue("\n \n Total number of species in data: {nrow(sp_sum)}"))

print(glue("\n \n Individuals per species:
           {sp_sum %>% print( n=36, max_footer_lines=0)}"))

print(glue("\n ##########################################################"))

print(glue("\n \n Total number of fixes in data: {sum(fix_sum$n)}"))

print(glue("\n \n Fixes per species:
           {evt_out %>% group_by(species) %>% summarize(n = n()) %>% as_tibble() %>% print( n=36, max_footer_lines=0)}"))

sink()

message("Script Complete!")
