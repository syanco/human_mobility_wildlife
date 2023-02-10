#!/usr/bin/env Rscript
#
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
  # rm(list=ls())
  library(here)
  
  .wd <- '/gpfs/ysm/project/jetz/ryo3/projects/covid'
  .test <- TRUE
  # rd <- here::here
  
  .datPF <- file.path(.wd,'raw_data/')
  .outPF <- file.path(.wd,"out/")
  
} else {
  library(docopt)
  library(rprojroot)
  
  .wd <- getwd()
  .script <-  thisfile()
  # rd <- is_rstudio_project$make_fix_file(.script)
  
  .datPF <- file.path(.wd,'raw_data/')
  .outPF <- file.path(.wd,"out/")
}

source(file.path(.wd,'analysis/src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(sf)
    library(raster)
    library(data.table)
  }))

# read in census block group geometries
message("reading in census block group geometries...")
cbg <- st_read(paste0(.datPF,"safegraph_open_census_data_2010_to_2019_geometry/cbg.geojson"))

ghm <- raster(paste0(.datPF,"gHM/gHM.tif"))

ghm <- raster::crop(ghm,cbg)

cbg_ghm <- st_transform(cbg, crs = st_crs(ghm))

cbg_ghm$ghm <- raster::extract(ghm,
                               cbg_ghm,
                               fun = mean,
                               na.rm = TRUE,
                               weights = TRUE)

st_write(cbg_ghm, paste0(.outPF,"safegraph_geometry_gHM/cbg_annotation.shp"))