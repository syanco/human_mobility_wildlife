#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# This script checks for data completeness and trims data absed on minimum data requirements


'
Filters data to minimums, reports out sample sizes

Usage:
filter_data_mins.r [--db=<db>] <wkmin> <minsp>
filter_data_mins.r (-h | --help)

Options:
-h --help     Show this screen.
-v --version     Show version.
-d --db=<db>  Path to movement database. Defaults to <wd>/data/move.db
-w --wkmin    Minimum sample size per week
-m --minsop    Minimum individuals per species
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
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

#---- Perform analysis ----#

#-- Load Cleaned Data

evt_cln <- tbl(db,'event_clean') %>%  collect()
ind_cln <- tbl(db, "individual_clean") %>%  collect()
# ind_cln <- tbl(db, "individual_trim") %>%  collect()
std_cln <- tbl(db, "study_clean") %>%  collect()
# std_cln <- tbl(db, "study_trim") %>%  collect()

#-- Remove species with too few individuals

# get ind count per species
sp_sum <- evt_cln %>%
  left_join(ind_cln %>% select(individual_id, taxon_canonical_name), by = "individual_id") %>% 
  mutate(species = taxon_canonical_name,
         ind_f = as.factor(individual_id)) %>% 
  group_by(species) %>%
  summarize(nind = length(unique(ind_f))) %>%
  filter(nind > .minsp) #require a minimum# of individuals

