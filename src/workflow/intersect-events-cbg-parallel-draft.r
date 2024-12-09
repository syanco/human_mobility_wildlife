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
  .datPF <- file.path(.wd,'raw_data/covid_movement_full_repo/raw_data/')
  .outPF <- file.path(.wd,'out/event-cbg-intersection/')
}

source(file.path(.wd,'src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(data.table)
    library(sf)
    library(foreach)
    library(doParallel)
    library(dplyr)
    library(assertthat)
  }))

# dir.create(.outPF, showWarnings = FALSE, recursive = TRUE)

# set up parallel backend
num_cores <- 5  # 40 cores per node on HPC "Pod"
cl <- makeCluster(num_cores)
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)

#---- Collect arguments ----#
# args = commandArgs(trailingOnly = TRUE)

# start_ix <- as.numeric(args[1])
# end_ix <- as.numeric(args[2])
# n <- as.numeric(args[3])

#---- Initialize database ----#

invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

# read in census block group geometries
message("reading in census block group geometries...")
cbg_sf <- read_sf(paste0(.datPF,"safegraph_open_census_data_2010_to_2019_geometry/cbg.geojson"))

gc()

# read in event table
message("reading in event table...")

evt_sf <- dbGetQuery(db,'SELECT event_id,lat,lon from event_final') %>%
  st_as_sf(coords = c("lon", "lat"), crs="+proj=longlat +datum=WGS84")

gc()

# calculate chunk sizes for parallelization
total_rows <- nrow(evt_sf)
chunk_size <- ceiling(total_rows / num_cores)
chunks <- split(1:total_rows, ceiling(seq_along(1:total_rows)/chunk_size))

message(sprintf("Processing %d rows in %d chunks using %d cores...", 
               total_rows, length(chunks), num_cores))

tmp_dir <- file.path(.outPF, "tmp_chunks")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# create logging file for parallel processing
progress_file <- file.path(.outPF, "progress.txt")
cat("Starting processing\n", file = progress_file)

message("intersecting events with census block groups...")

#evt_sf <- dbGetQuery(db,'SELECT * from event_clean') %>%
#  collect() %>%
#  st_as_sf(coords = c("lon", "lat"), crs="+proj=longlat +datum=WGS84")

# no need for the following line, since not using DSQ package
# evt_sf <- evt_sf[start_ix:end_ix,]

tryCatch({
  results <- foreach(
    chunk = chunks,
    chunk_id = seq_along(chunks),
    .packages = c('sf', 'dplyr', 'data.table'),
    .errorhandling = 'pass'  # continue processing if one chunk fails
  ) %dopar% {

    chunk_msg <- sprintf("Starting chunk %d of %d\n", chunk_id, length(chunks))
    cat(chunk_msg, file=progress_file, append=TRUE)

    evt_chunk <- evt_sf[chunk,]

    # intersect event table with census block group geometries
    evt_cbg <- st_intersection(evt_chunk,cbg_sf) %>%
      rename(cbg_2010 = CensusBlockGroup) %>%
      st_drop_geometry()

    gc()

    # write out new table with annotations
    chunk_file <- file.path(tmp_dir, sprintf("chunk_%02d.csv", chunk_id))
    # TODO: add step where it checks if the out/even-cbg-intersection subdir exists, and creates it first if not 
    fwrite(evt_cbg, chunk_file)

    rm(evt_cbg)
    gc()

    chunk_msg <- sprintf("Completed chunk %d of %d\n", chunk_id, length(chunks))
    cat(chunk_msg, file=progress_file, append=TRUE)

    # return NULL since we're writing files
    NULL
  }
}, error = function(e) {

  cat(sprintf("Error occurred: %s\n", e$message), file=progress_file, append=TRUE)
  message("Error; Check progress.txt for details.")
  
}, finally = {
  # Don't delete tmp_dir in case we need to recover data
  cat("Processing finished or interrupted. Temporary files preserved in case needed.\n", 
      file=progress_file, append=TRUE)
})


stopCluster(cl)

# concatenate chunk files into one CSV
message("Concatenating chunk files...")
output_file <- file.path(.outPF, "event-cbg-intersection.csv")

# write headers from first file
system(sprintf("head -n 1 %s/chunk_01.csv > %s", temp_dir, output_file))
# append data from all files (skipping headers)
system(sprintf("for f in %s/chunk_*.csv; do tail -n +2 \"$f\"; done >> %s", 
              tmp_dir, output_file))

dbDisconnect(db)

message("done!")