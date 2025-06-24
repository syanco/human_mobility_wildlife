#!/usr/bin/env Rscript
# DESCRIPTION #
#
# This script determines the level of global human mofication (GHM) for animal 
# movement locations for the COVID-19 Animal Movement Project. GHM data is sourced
# from an existing global map in raster format, which estimates cumulative landscape 
# modification within 1-km2 pixels from 13 human stressor datasets.
#
# See project documentation for details about anticipated directory structure.

# ==== Breezy setup ====

#---- Input Parameters ----#
if(interactive()) {

  library(here)
  
  .wd <- '/home/sy522/project/covid-19_movement'
  .test <- TRUE
  .dbPF <- '/home/sy522/project/covid-19_movement/processed_data/mosey_mod.db'
  .datPF <- file.path(.wd,'raw_data/')
  .outPF <- file.path(.wd,'out/')
  
} else {

  library(docopt)
  library(rprojroot)

  .wd <- getwd()
  .dbPF <- '/tmp/mosey_mod.db'
  .datPF <- file.path('/home/julietcohen/covid_movement_full_repo/raw_data/')
  .outPF <- file.path(.wd,'out/')
}

source(file.path(.wd,'src/startup.r'))

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
evt_sf <- dbGetQuery(db,'SELECT event_id,lon,lat from event_final') %>%
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

# write out new table with annotations
message("writing out new event table...")
fwrite(evt_ghm, paste0(.outPF, "event-annotation/event_ghm.csv"))

dbDisconnect(db)

message("intersection with human modification done!")