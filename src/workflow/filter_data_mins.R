#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# This script checks for data completeness and trims data absed on minimum data requirements


'
Check for data completeness and trims data absed on minimum data requirements

Usage:
filter_data_mins.r [--db=<db>] 
filter_data_mins.r (-h | --help)

Options:
-h --help     Show this screen.
-v --version     Show version.
-d --db=<db> Path to movement database. Defaults to <wd>/data/move.db
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/covid-19_movement'
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



# dbBegin(db)

#---- Perform analysis ----#

#-- Clean Outliers

evt_trm <- tbl(db,'event_clean') %>%  collect()
ind_trm <- tbl(db, "individual_clean") %>%  collect()
std_trm <- tbl(db, "study_clean") %>%  collect()