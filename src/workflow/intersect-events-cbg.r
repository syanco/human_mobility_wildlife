#!/usr/bin/env Rscript --vanilla
#
# DESCRIPTION #
#
# This script determines the administrative units for records for the COVID-19 Animal Movement Project
# See project documentation for details about anticipated directory structure.
#
# Major tasks fof this script:
#   * Annotate event dataset with:
#     * CensusBlockGroup (12 digit FIPS code)
#     * BlockGroup (1 digit FIPS code)
#     * TractCode (6 digit FIPS)
#     * CountyFIPS (3 digit FIPS)
#     * StateFIPS (2 digit FIPS)
#     * County (character string)
#     * State (2 character code)
#   *write out csv with administrative info
#
# This script implements the breezy philosophy: github.com/benscarlson/breezy


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
  
  .wd <- getwd()
  .test <- TRUE
  rd <- here::here
  
  .dbPF <- file.path(.wd, 'processed_data/mosey_mod.db')
  # .dbPF <- '/home/sy522/project/covid-19_movement/processed_data/mosey_swap_mod.db'
  .datPF <- file.path(.wd,'raw_data/')
  .outPF <- file.path(.wd,'analysis/event-cbg-intersection/')
  
} else {
  library(docopt)
  library(rprojroot)
  
  .wd <- getwd()
  # .script <-  thisfile()
  # rd <- is_rstudio_project$make_fix_file(.script)
  
  .dbPF <- file.path('/tmp/mosey_mod.db')
  # .dbPF <- '/home/sy522/project/covid-19_movement/processed_data/mosey_swap_mod.db'
  # .datPF <- file.path(.wd,'raw_data/covid_movement_full_repo/raw_data/')
  .datPF <- file.path('/home/julietcohen/covid_movement_full_repo/raw_data/')
  .outPF <- file.path(.wd,'out/event-cbg-intersection/')
}

source(file.path(.wd,'src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(data.table)
    library(sf)
    library(assertthat)
    library(dplyr)
  }))

#---- Initialize database ----#

invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

# read in census block group geometries
message("reading in census block group geometries...")
cbg_sf <- read_sf(paste0(.datPF,"safegraph_open_census_data_2010_to_2019_geometry/cbg.geojson"))

# read in event table
message("reading in event table...")

evt_sf <- dbGetQuery(db,'SELECT event_id,lat,lon from event_final') %>%
  st_as_sf(coords = c("lon", "lat"), crs="+proj=longlat +datum=WGS84")

#evt_sf <- dbGetQuery(db,'SELECT * from event_clean') %>%
#  collect() %>%
#  st_as_sf(coords = c("lon", "lat"), crs="+proj=longlat +datum=WGS84")

# intersect event table with census block group geometries
message("intersecting events with census block groups...")
evt_cbg <- st_join(evt_sf, 
                   cbg_sf, 
                   join = st_intersects, 
                   left = FALSE) %>% # do not retain points that do not intersect any CBG
          rename(cbg_2010 = CensusBlockGroup) %>%
          st_drop_geometry() %>%
          # if a point falls within multiple cbg polygons, keep the first observation
          group_by(event_id) %>%
          slice(1) %>%
          ungroup()

# write out new table with annotations
# TODO: add step where it checks if the out/even-cbg-intersection subdir exists, and creates it first if not 
message("writing out csv...")
fwrite(evt_cbg, paste0(.outPF,"event-cbg-intersection.csv"))

dbDisconnect(db)

message("done!")