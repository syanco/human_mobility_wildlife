#!/usr/bin/env Rscript 
#
# DESCRIPTION #
#
# This script links human mobility data with animal events for the COVID-19 Animal Movement Project
# Human mobility sourced from SafeGraph, provided by Song Gao
#
# See project documentation for details about anticipated directory structure.
# This script implements the breezy philosophy: github.com/benscarlson/breezy
#
# Major tasks fof this script:
#   * Annotate event dataset with:
#     * daily and hourly device counts
#   *write out table with decive counts

# ==== Breezy setup ====

#---- Input Parameters ----#
if(interactive()) {

  library(here)
  
  .wd <- '/gpfs/ysm/project/jetz/ryo3/projects/covid'
  .test <- TRUE
  .dbPF <- '/home/sy522/project/covid-19_movement/processed_data/mosey_mod.db'
  .datPF <- file.path(.wd,'out/')
  .pdPF  <- file.path(.wd,'processed_data')
  .outPF <- file.path(.wd,'out/')
  
} else {
  
  .wd <- getwd()
  .dbPF <- '/tmp/mosey_mod.db'
  .datPF <- file.path(.wd,'out/')
  .pdPF  <- file.path('/home/julietcohen/covid_movement_full_repo/processed_data/')
  .outPF <- file.path(.wd,'out/')
}

message("start safegraph annotation")
source(file.path(.wd,'src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(data.table)
  }))

#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

message("reading in event table...")
evt_df <- dbGetQuery(db,'SELECT event_id,timestamp from event_final') %>%
  separate(timestamp, c("date",NA), sep = " ", remove = FALSE) %>%
  mutate("date_hour" = str_trunc(timestamp,13,"right","")) %>%
  collect()

evt_cbg <- fread(paste0(.datPF, "event-annotation/event_cbg.csv"), colClasses = "character")

reformatted_files_daily <- list.files(paste0(.pdPF,"safegraph/counties-dates-2-10-22-reformatted/daily-data"), full.names = TRUE)

# combine all data
message("reading in safegraph data...")
daily_data <- data.table::rbindlist(lapply(reformatted_files_daily, data.table::fread, colClasses = "character"),use.names = TRUE) %>%
  select(cbg,date,count) %>%
  rename(safegraph_daily_count = count) %>%
  mutate(cbg_2010 = as.character(cbg),
         date = as.character(date)) %>%
  select(-cbg)

message("joining events with safegraph data...")
evt_sg <- evt_df %>%
  mutate(event_id = as.character(event_id)) %>%
  left_join(.,evt_cbg, by = "event_id") %>%
  left_join(.,daily_data, by = c("cbg_2010", "date")) %>%
  select(-date)

message("writing out new event table...")
fwrite(evt_sg, paste0(.outPF, "event-annotation/event_sg.csv"))

dbDisconnect(db)

message("safegraph annotation done!")