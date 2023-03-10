#!/usr/bin/env Rscript 
#
#CONDA: covid
#
# This scriptMerges the new and replacement studies in mosey_swap_mod.db into the 
# main analysis db (mosey_mod.db)

'
Merge mosey_swap_mod.db with mosey_mod.db

Usage:
merge_dbs.r 
merge_dbs.r (-h | --help)

Parameters:
  db_swap: path to the swap databse. 
  db_main: path to the main database.

Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/covid-19_movement'
  rd <- here::here
  
  # .dbmain <- '~/project/covid-19_movement/processed_data/mosey_mod.db'
  .dbmain <- '~/projects/covid-19_movement/processed_data/mosey_mod_20220303.db'
  .dbswap <- '~/projects/covid-19_movement/processed_data/mosey_swap_mod.db'
  
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  # .dbmain <- '~/project/covid-19_movement/processed_data/mosey_mod.db'
  .dbmain <- '~/project/covid-19_movement/processed_data/mosey_mod_20220303.db'
  .dbswap <- '~/project/covid-19_movement/processed_data/mosey_swap_mod.db'

  
}

#---- Initialize Environment ----#
# Start time
t0 <- Sys.time()

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
    library(tidyverse)
  }))

# Load Control Files
# ctfnew <- read_csv('~/project/covid-19_movement/ctfs/new_studies.csv')
# ctfswap <- read_csv('~/project/covid-19_movement/ctfs/swaps.csv')
ctfnew <- read_csv('~/projects/covid-19_movement/ctfs/new_studies.csv')
ctfswap <- read_csv('~/projects/covid-19_movement/ctfs/swaps.csv')

swapstudies <- ctfswap$Old_study_name # vec of studies to swap out (not including the new ones)

# Run startup
source(file.path(.wd,'analysis/src/startup.r'))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbmain)))
dbmain <- dbConnect(RSQLite::SQLite(), .dbmain, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(dbmain))>0))

invisible(assert_that(file.exists(.dbswap)))
dbswap <- dbConnect(RSQLite::SQLite(), .dbswap, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(dbswap))>0))


# merge individual tables
message("Merging individual tables...")
ind_swap <- tbl(dbswap, 'individual') %>% 
  collect()

ind_main <- tbl(dbmain,'individual') %>% # load main db event table
  collect() %>% # collect into memory  
  filter(study_id %notin% swapstudies) %>% # remove studies that were swapped
  bind_rows(ind_swap)   # bind new event table in

swapinds <- tbl(dbmain, 'individual') %>% 
  filter(study_id %in% swapstudies) %>% 
  collect() %>% 
  pull(individual_id)

# write back to db
message("Writing individual table back to db...")
dbWriteTable(dbmain, 'individual', ind_main, overwrite = T)

# merge event tables
message("Merging event tables...")
evt_swap <- tbl(dbswap, 'event_clean') %>% 
  collect()

evt_main <- tbl(dbmain,'event_clean') %>% # load main db event table
  collect() %>% # collect into memory  
  filter(individual_id %notin% swapinds) %>% # remove studies that were swapped
  bind_rows(evt_swap)   # bind new event table in

  # write back to db
message("Writing event table back to db...")
dbWriteTable(dbmain, 'event_clean', evt_main, overwrite = T)

# merge study tables
message("Merging study tables...")

std_swap <- tbl(dbswap, 'study') %>% 
  collect()

std_main <- tbl(dbmain,'study') %>% # load main db event table
  collect() %>% # collect into memory  
  filter(study_id %notin% swapstudies) %>% # remove studies that were swapped
  bind_rows(std_swap)   # bind new event table in

# write back to db
message("Writing study table back to db...")
dbWriteTable(dbmain, 'study', std_main, overwrite = T)

message("Disconnecting from db...")
dbDisconnect(dbmain)
dbDisconnect(dbswap)

# Export tables as csvs
# message("Exporting data to csv...")
# write_csv(evt_main, file = "out/data_export/events.csv")
# write_csv(std_main, file = "out/data_export/study.csv")
# write_csv(ind_main, file = "out/data_export/individual.csv")


#---- Finalize script ----#

message(glue('Script complete in {diffmin(t0)} minutes'))
