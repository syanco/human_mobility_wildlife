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
  # rm(list=ls())
  library(here)
  
  .wd <- '/gpfs/ysm/project/jetz/ryo3/projects/covid'
  .test <- TRUE
  .dbPF <- '/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db'
  .datPF <- file.path(.wd,'raw_data/')
  .outPF <- file.path(.wd,"out/")
  
} else {
  library(docopt)
  library(rprojroot)
  
  .wd <- getwd()
  .script <-  thisfile()
  .dbPF <- '/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db'
  .datPF <- file.path(.wd,'raw_data/')
  .outPF <- file.path(.wd,"out/")
}

message("start census annotation")

source(file.path(.wd,'analysis/src/startup.r'))

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
evt_df <- dbGetQuery(db,'SELECT * from event_cbg') %>%
  collect()

message("reading in census data...")

acs2019 <- fread(paste0(.datPF,"safegraph_open_census_data_2019/data/cbg_b01.csv")) %>%
  select(census_block_group,B01003e1) %>%
  rename(cbg_2010 = census_block_group,
         total_population_2019 = B01003e1) %>%
  mutate(cbg_2010 = as.character(cbg_2010))

message("joining events with census data...")
evt_census <- left_join(evt_df,acs2019, by = c("cbg_2010" = "cbg_2010"))

message("writing out new event table...")
fwrite(evt_census, paste0(.outPF, "event-annotation/event_census.csv"))

dbDisconnect(db)


message("census annotation done!")