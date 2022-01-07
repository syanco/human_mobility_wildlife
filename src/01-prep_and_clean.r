#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# This script preps and cleans data for the COVID-19 Animal Movement Project
# See project documentation for details about anticipated directory structure.
# Major tasks of this script:
#   * Annotate dataset with the study periods:
#     * Pre-LD 2019
#     * LD 2019
#     * Post-LD 2019
#     * Pre-LD 2020
#     * LD 2020
#     * Post-LD 2020
#   * Remove events outside relevant study periods.
#   * Basic data cleaning:
#     * ID and remove outliers


'
Usage:
01-prep_and_clean.r [--db=<db>] [--cores=<cores>]
01-prep_and_clean.r (-h | --help)

Control files:
  src/ctfs/dates.csv


Options:
-h --help     Show this screen.
-v --version     Show version.
-d --db=<db> Path to movement database. Defaults to <wd>/data/move.db
-c --cores=<cores>  The number of cores
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/covid-19_movement'
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
  .nc <- 2
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  .dbPF <- makePath(ag$db)
  .nc <- ag$cores

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
    library(atlastools)
    library(foreach)
    library(doMC)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)


#---- Load control files ----#
periods <- read_csv(file.path(.wd,'analysis/ctfs/dates.csv'),
                    col_types=list("date" = col_date(format = "%m/%d/%Y"))) 

#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

evt <- tbl(db,'event')
# dbBegin(db)
#---- Perform analysis ----#

#-- Make a filtered table by study period

# Code below can be used to crib the SQL query from `dbplyr`...

# evt %>% 
#   filter((timestamp > t1 & timestamp < t2) | (timestamp > t3 & timestamp < t4)) %>% 
#   mutate(yr = strftime('%Y', timestamp),
#                        trt = case_when(
#                          timestamp >= !!periods$date[periods$cutpoint == "start_pre-ld_2019"] &
#                            timestamp < !!periods$date[periods$cutpoint == "start_ld_2019"] ~ "pre-ld_2019",
#                          timestamp >= !!periods$date[periods$cutpoint == "start_ld_2019"] &
#                            timestamp < !!periods$date[periods$cutpoint == "start_post-ld_2019"] ~ "ld_2019",
#                          timestamp >= !!periods$date[periods$cutpoint == "start_post-ld_2019"] &
#                            timestamp < !!periods$date[periods$cutpoint == "stop_2019"] ~ "post-ld_2019",
#                          timestamp >= !!periods$date[periods$cutpoint == "start_pre-ld_2020"] &
#                            timestamp < !!periods$date[periods$cutpoint == "start_ld_2020"] ~ "pre-ld_2020",
#                          timestamp >= !!periods$date[periods$cutpoint == "start_ld_2020"] &
#                            timestamp < !!periods$date[periods$cutpoint == "start_post-ld_2020"] ~ "ld_2020",
#                          timestamp >= !!periods$date[periods$cutpoint == "start_post-ld_2020"] &
#                            timestamp < !!periods$date[periods$cutpoint == "stop_2020"] ~ "post-ld_2020")
#                        ) %>% show_query()

# SQLite query to filter out events outside study periods, add a year variable,
# and a study period variable.  Uses ctf to assign cutpoints
q <- glue("
  CREATE TABLE event_mod AS
    SELECT *, strftime('%Y', `timestamp`) AS `yr`, 
      CASE
        WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_pre-ld_2019']}' 
          AND `timestamp` < '{periods$date[periods$cutpoint == 'start_ld_2019']}') 
            THEN ('pre-ld_2019')
        WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_ld_2019']}' 
          AND `timestamp` < '{periods$date[periods$cutpoint == 'start_post-ld_2019']}') 
            THEN ('ld_2019')
        WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_post-ld_2019']}' 
          AND `timestamp` < '{periods$date[periods$cutpoint == 'stop_2019']}') 
            THEN ('post-ld_2019')
        WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_pre-ld_2020']}' 
          AND `timestamp` < '{periods$date[periods$cutpoint == 'start_ld_2020']}') 
            THEN ('pre-ld_2020')
        WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_ld_2020']}' 
          AND `timestamp` < '{periods$date[periods$cutpoint == 'start_post-ld_2020']}') 
            THEN ('ld_2020')
        WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_post-ld_2020']}' 
          AND `timestamp` < '{periods$date[periods$cutpoint == 'stop_2020']}') 
            THEN ('post-ld_2020')
      END AS `trt`
    FROM `event`
      WHERE ((`timestamp` > '{periods$date[periods$cutpoint == 'start_pre-ld_2019']}' 
        AND `timestamp` < '{periods$date[periods$cutpoint == 'stop_2019']}') 
        OR (`timestamp` > '{periods$date[periods$cutpoint == 'start_pre-ld_2020']}' 
        AND `timestamp` < '{periods$date[periods$cutpoint == 'stop_2020']}'))"
  )

# Execute the query
dbExecute(db, q)

# dbRollback(db)
# dbCommit(db)


#-- Clean Outliers

# Get list of inds
inds <- tbl(db, 'event_mod') %>% 
  select(individual_id) %>% 
  distinct() %>% 
  collect() %>% 
  c()

registerDoMC(.nc)
foreach(i = 1:length(inds)) %dopar% {
  
  #TODO: fill in with atlas_tools (page 9)
  
  )

#---- Finalize script ----#

dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))