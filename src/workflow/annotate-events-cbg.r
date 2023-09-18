#!/usr/bin/env Rscript --vanilla
#
# DESCRIPTION #
#
# This script links census geometry info with animal events for the COVID-19 Animal Movement Project
# Event table is spatially intersected with census block group (cbg) geometries in intersect-events-cbg.R
# Area for each geometry is computed in compute-cbg-area.R
#
# See project documentation for details about anticipated directory structure.
# This script implements the breezy philosophy: github.com/benscarlson/breezy
#
# Major tasks of this script:
#   * combine event/cbg intersection files
#   * read in cbg area
#   * combine event table with associated cbgs and cbg area

# ==== Breezy setup ====

#'
#Template
#Usage:
#script_template <taxa> <dat> <out> 
#script_template (-h | --help)
#Parameters:
#  dat: path to input csv file. 
#  out: path to output directory.
#Options:
#-h --help     Show this screen.
#-v --version     Show version.
#' -> doc

#---- Input Parameters ----#
if(interactive()) {
  rm(list=ls())
  library(here)
  
  .wd <- '/gpfs/ysm/project/jetz/ryo3/projects/covid'
  .test <- TRUE
  # rd <- here::here
  
  .dbPF <- '/home/sy522/project/covid-19_movement/processed_data/mosey_mod.db'
  # .dbPF <- '/home/sy522/project/covid-19_movement/processed_data/mosey_swap_mod.db'
  .datPF <- file.path(.wd,'out/')
  .outPF <- file.path(.wd,'out/')
  
} else {
  # library(docopt)
  # library(rprojroot)
  
  .wd <- getwd()
  # .script <-  thisfile()
  # rd <- is_rstudio_project$make_fix_file(.script)
  .dbPF <- '/home/sy522/project/covid-19_movement/processed_data/mosey_mod_2023.db'
  # .dbPF <- '/home/sy522/project/covid-19_movement/processed_data/mosey_swap_mod.db'
  .datPF <- file.path(.wd,'out/')
  .outPF <- file.path(.wd,'out/')
}

message("start safegraph annotation")
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

message("reading in files...")
files <- list.files(paste0(.datPF,"event-cbg-intersection/"),pattern = "*.csv",full.names = TRUE)
intersection = data.table::rbindlist(lapply(files, data.table::fread, colClasses = "character"),use.names = TRUE)

data_counties <- intersection %>%
  unite(col = "county", StateFIPS:CountyFIPS, sep ="") %>%
  distinct(county)

data_cbg <- intersection %>%
  distinct(cbg_2010) %>%
  rename(CensusBlockGroup = cbg_2010)

fwrite(data_counties, paste0(.outPF,"safegraph-summary/counties-list.csv"))
fwrite(data_cbg, paste0(.outPF,"safegraph-summary/census-block-group-list.csv"))

area <- fread(paste0(.datPF,"event-annotation/cbg-area.csv"), colClasses = "character") %>%
  select(cbg_2010, cbg_area_m2)

message("reading in event table...")
evt_df <- dbGetQuery(db,'SELECT event_id from event_final')

message("joining event table with cbg info..")
evt_cbg <- evt_df %>%
  mutate(event_id = as.character(event_id)) %>%
  left_join(., intersection, by = "event_id") %>%
  left_join(., area, by = "cbg_2010")

message("writing out new event table...")
fwrite(evt_cbg, paste0(.outPF, "event-annotation/event_cbg.csv"))
#dbWriteTable(conn = db, name = "event_cbg", value = evt_cbg, append = FALSE, overwrite = T)

dbDisconnect(db)

message("cbg annotation done!")