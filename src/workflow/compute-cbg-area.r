#!/usr/bin/env Rscript --vanilla
#
# DESCRIPTION #
#
# This script determines the administrative units for records for the COVID-19 Animal Movement Project
# See project documentation for details about anticipated directory structure.
#
# Major tasks of this script:
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
  
  .dbPF <- '/home/sy522/project/covid-19_movement/processed_data/mosey_mod.db'
  # .dbPF <- '/home/sy522/project/covid-19_movement/processed_data/mosey_swap_mod.db'
  .test <- TRUE
  # rd <- here::here
  
  .datPF <- file.path(.wd,'data/')
  .outPF <- file.path(.wd,'out/event-annotation/')
  
} else {
  library(docopt)
  library(rprojroot)
  
  .wd <- getwd()
  # .script <-  thisfile()
  # rd <- is_rstudio_project$make_fix_file(.script)
  
  .datPF <- file.path(.wd,'raw_data/')
  .outPF <- file.path(.wd,'out/event-annotation/')
}

source(file.path(.wd,'analysis/src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(data.table)
    library(sf)
  }))


#---- start analysis ----#

# read in census block group geometries
message("reading in census block group geometries...")
cbg_sf <- st_read(paste0(.datPF,"safegraph_open_census_data_2010_to_2019_geometry/cbg.geojson"))

message("calculating area...")
cbg_sf$cbg_area_m2 <- st_area(cbg_sf)

cbg <- cbg_sf %>%
  rename(cbg_2010 = CensusBlockGroup) %>%
  st_drop_geometry()

message("writing out results...")
fwrite(cbg, paste0(.outPF,"cbg-area.csv"))

message("done!")

