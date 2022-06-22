#!/usr/bin/env Rscript #
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
#   *write out table with administrative info
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
  
  .wd <- '/gpfs/ysm/project/jetz/ryo3/projects/covid'
  .test <- TRUE
  rd <- here::here
  
  .dbPF <- '/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db'
  .datPF <- file.path(.wd,'data/')
  .outPF <- file.path(.wd,'analysis/')
  
} else {
  library(docopt)
  library(rprojroot)
  
  .wd <- '/gpfs/ysm/project/jetz/ryo3/projects/covid'
  .script <-  thisfile()
  rd <- is_rstudio_project$make_fix_file(.script)
  
  .dbPF <- '/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db'
  .datPF <- file.path(.wd,'data/')
  .outPF <- file.path(.wd,'analysis/')
}

source(file.path(.wd,'/src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(data.table)
    library(raster)
    library(sf)
  }))

#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

# read in global human modification layer
message("reading in human modification...")
ghm <- raster(paste0(.datPF,"gHM/gHM.tif"))

# read in event table
message("reading in event table...")
evt_sf <- dbGetQuery(db,'SELECT event_id,lon,lat from event_clean') %>%
  collect() %>%
  st_as_sf(coords = c("lon", "lat"), crs="+proj=longlat +datum=WGS84")

# transform to raster reference system
message("transform event table...")
evt_sf <- st_transform(evt_sf,st_crs(ghm))

# extract values at points
message("intersecting events with human modification...")
evt_sf$ghm <- raster::extract(ghm,evt_sf)

evt_ghm <- evt_sf %>%
  st_drop_geometry()

head(evt_ghm)

# write out new table with annotations
message("writing out new event table...")
fwrite(evt_ghm, paste0(.outPF, "event-annotations/event_ghm.csv"))
#dbWriteTable(conn = db, name = "event_ghm", value = evt_ghm, append = FALSE, overwrite = T)

dbDisconnect(db)

message("intersection with human modification done!")