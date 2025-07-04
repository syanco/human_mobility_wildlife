#!/usr/bin/env Rscript 
#
#CONDA: covid
#
# This script uses previously calculated dBBMMs to assess animal space use 
# during COVID-19 lockdowns

# TODO: Clean the script - unused object created, poorly commented, etc
# 
# ==== Setup ====

'
Calculate space use before and during COVID-19 lockdowns using previously 
estimated dBBMMs

Usage:
calc-space-use.r <out> <db> <ctf> <nc> [-c]
calc-space-use.r (-h | --help)

Parameters:
  --out: path to directory storing dBBMM outputs. 
  --db: path to dtabase
  --ctf: path to dbbmm log file
  --nc: number of cores for parallel processing

  
Options:
-h --help     Show this screen.
-v --version     Show version.
-c --continue   Indicates if script should first check for previous entries before calculating metrics
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '/home/julietcohen/repositories/human_mobility_wildlife'
  .outPF <- file.path(.wd,'out_size_test')
  .dbPF <- file.path(.wd,'processed_data/intermediate_db_copies/mosey_mod_clean-movement_complete.db')
  .ctf <- file.path(.wd, "out/dbbmm_log.csv")
  .continue = T
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'src/funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  .ctf <- makePath(ag$ctf)
  .nc <- ag$nc
  .continue <- as.logical(ag$continue)
  
}

#---- Initialize Environment ----#
# Start time
t0 <- Sys.time()

# Run startup
source(file.path(.wd,'src/startup.r'))

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
    library(raster)
    library(move)
    library(doMC)
    library(foreach)
    library(sf)
    library(tidyverse )
    library(glue)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'), full.names=TRUE) %>%
  walk(source)


#---- Initialize database ----#
message("Initializing database connection and control files...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

ind <- tbl(db,'individual_clean') %>% 
  collect() %>% 
  pull(individual_id)

evt_anno <- tbl(db, "event_clean") %>% 
  collect()

dbDisconnect(db)

# yearvec <- c("2019", "2020")
# trtvec <- c("pre-ld", "ld")

# SG from csv
evt_sg <- read.csv("out/event-annotation/event_sg.csv", 
                   colClasses = c(event_id = "integer64"))

# GHM from csv
evt_ghm <- read.csv("out/event-annotation/event_ghm.csv", 
                    colClasses = c(event_id = "integer64"))

ctf <- read.csv(.ctf)
#TODO: rm below after a real run with log
ctf <- ctf[!duplicated(ctf),] %>% 
  filter(produced == 1)

# #- get the phylo control covar
# #TODO
# sp2019 <- raster("raw_data/geoserver_1643663843611/data.tif")
# sp2020 <- raster("raw_data/geoserver_1643663856740/data.tif")

#---- Perform analysis ----#
message("Gathering movement data...")

# Load existing file in cases where you want to pick upn after eg a timeout
existing_out <- read_csv(glue("{.outPF}/dbbmm_size.csv"))
# existing_out <- existing_out[0,]
#+++++++++++++++++++++#
message(glue("Starting {.nc} cores"))
registerDoMC(.nc)

# Toggle `%do%` to `%dopar%` for HPC, %do% for local
foreach(i = 1:nrow(ctf), .errorhandling = "pass", .inorder = F) %dopar% {
  
  out_check <- existing_out %>% 
    filter(ind_id == glue("{ctf$ind_id[i]}"))
  
  if(.continue){
    if(nrow(out_check) == 0){
      message(glue("Starting {ctf$species[i]} ind {ctf$ind_id[i]}, year {ctf$year[i]}"))
      tryCatch({
        load(glue("{.outPF}/dbbmms/dbbmm_{ctf$ind_id[i]}_{ctf$year[i]}.rdata"))
        
        r <- tmp_out$`dBBMM Object`
        # r
        # plot(sqrt(r))
        rb <- UDStack(r)
        UDr <- getVolumeUD(rb)
        # plot(UDr)
        for(j in 1:nlayers(UDr)){
          
          # Get UD area:
          # 1. create a binary raster with threshold at 95%
          ud95 <- UDr[[j]]<=.95
          # 2. calculate total area in square meters (UD units are already meters)
          a <- sum(values(ud95))*res(ud95)[1]*res(ud95)[1]
          
          # get week number
          week <- as.numeric(substring(names(UDr[[j]]),2))
          
          # Get event IDs underlying dBBMM (used as key to filter annotations below)
          evtids <- tmp_out$events %>%
            filter(wk == week) %>% 
            pull(event_id)
          
          # Get safegraph data
          sg <- evt_sg %>% 
            filter(event_id %in% evtids) %>% 
            summarize(sg = mean(safegraph_daily_count, na.rm = T))
          
          # Get area of CBG to normalize safegraph data
          cbg_area <- evt_sg %>% 
            filter(event_id %in% evtids) %>% 
            summarize(cbg = mean(cbg_area_m2, na.rm = T))
          
          # Get GHM data
          ghm <- evt_ghm %>% 
            filter(event_id %in% evtids) %>% 
            summarize(ghm = mean(ghm, na.rm = T))
          
          # Get NDVI
          ndvi <- evt_anno %>% 
            filter(event_id %in% evtids) %>% 
            summarize(ndvi = mean(ndvi, na.rm = T))
          
          # Get TMAX
          tmax <- evt_anno %>% 
            filter(event_id %in% evtids) %>% 
            summarize(tmax = mean(tmax, na.rm = T))
          
          #unpack underlying data
          evt_tmp <- tmp_out$events 
          evt_mv <- move(x=evt_tmp$lon, y=evt_tmp$lat, time=evt_tmp$timestamp, 
                         trt = evt_tmp$trt,
                         proj=CRS("+proj=longlat"))
          evt_sf <- st_as_sf(evt_tmp, coords = c("lon", "lat"), crs = 4326)
          
          # get n pts
          n <- nrow(evt_tmp)
          
          # get area of bounding box of uinderlying points      
          a_bb <- st_area(st_make_grid(evt_sf, n=1))
          
          # get median fix rate
          fixmed <- median(timeLag(x=evt_mv, units="mins"))
          
          # get meann horiz accuracy
          m_error <- mean(na.omit(evt_tmp$horizontal_accuracy))
          
          # Write Out Results
          out <- matrix(c(ctf$species[i], 
                          ctf$ind_id[i], 
                          ctf$study_id[i], 
                          ctf$year[i], 
                          week, 
                          a, 
                          sg, 
                          ghm,
                          cbg_area, 
                          ndvi, 
                          tmax,
                          n, 
                          a_bb, 
                          fixmed, 
                          m_error),
                        nrow = 1)
          
          message(glue("Writing info for {ctf$species[i]} ind {ctf$ind_id[i]}, year {ctf$year[i]}, 
                   week {week}"))
          write.table(out, glue("{.outPF}/dbbmm_size.csv"), append = T, 
                      row.names = F, col.names = F, sep = ",")
          
        } #j
      }, error = function(e){cat(glue("ERROR: Size calulation failed for individual 
                                  {ctf$species[i]} {ctf$ind_id[i]}, year {ctf$year[i]}", 
                                  "\n"))})
    # if file hasn't been written
      } else {message(glue("Metrics for individual {ctf$ind_id[i]} already calculated and continue is set to T, gotta keep movin' on..."))}
  } else { # if continue is set to false, just do the thing
    # TODO:  this is a lot of duplicated code, should probably just make it a function but lazy...
    message(glue("Starting {ctf$species[i]} ind {ctf$ind_id[i]}, year {ctf$year[i]}"))
    tryCatch({
      load(glue("{.outPF}/dbbmms/dbbmm_{ctf$ind_id[i]}_{ctf$year[i]}.rdata"))
      
      r <- tmp_out$`dBBMM Object`
      # r
      # plot(sqrt(r))
      rb <- UDStack(r)
      UDr <- getVolumeUD(rb)
      # plot(UDr)
      for(j in 1:nlayers(UDr)){
        
        # Get UD area
        ud95 <- UDr[[j]]<=.95
        a <- sum(values(ud95))*res(ud95)[1]*res(ud95)[1]
        
        # get week number
        week <- as.numeric(substring(names(UDr[[j]]),2))
        
        # Get event IDs underlying dBBMM (used as key to filter annotations below)
        evtids <- tmp_out$events %>%
          filter(wk == week) %>% 
          pull(event_id)
        
        # Get safegraph data
        sg <- evt_sg %>% 
          filter(event_id %in% evtids) %>% 
          summarize(sg = mean(safegraph_daily_count, na.rm = T))
        
        # Get area of CBG to normalize safegraph data
        cbg_area <- evt_sg %>% 
          filter(event_id %in% evtids) %>% 
          summarize(sg = mean(cbg_area_m2, na.rm = T))
        
        # Get GHM data
        ghm <- evt_ghm %>% 
          filter(event_id %in% evtids) %>% 
          summarize(ghm = mean(ghm, na.rm = T))
        
        # Get NDVI
        ndvi <- evt_anno %>% 
          filter(event_id %in% evtids) %>% 
          summarize(ndvi = mean(ndvi, na.rm = T))
        
        # Get tmax
        tmax <- evt_anno %>% 
          filter(event_id %in% evtids) %>% 
          summarize(tmax = mean(tmax, na.rm = T))
        
        #unpack underlying data
        evt_tmp <- tmp_out$events 
        evt_mv <- move(x=evt_tmp$lon, y=evt_tmp$lat, time=evt_tmp$timestamp, 
                       trt = evt_tmp$trt,
                       proj=CRS("+proj=longlat"))
        evt_sf <- st_as_sf(evt_tmp, coords = c("lon", "lat"), crs = 4326)
        
        # get n pts
        n <- nrow(evt_tmp)
        
        # get area of bounding box of uinderlying points      
        a_bb <- st_area(st_make_grid(evt_sf, n=1))
        
        # get median fix rate
        fixmed <- median(timeLag(x=evt_mv, units="mins"))
        
        # get meann horiz accuracy
        m_error <- mean(na.omit(evt_tmp$horizontal_accuracy))
        
        # Write Out Results
        out <- matrix(c(ctf$species[i], 
                        ctf$ind_id[i], 
                        ctf$study_id[i], 
                        ctf$year[i], 
                        week, 
                        a, 
                        sg, 
                        ghm,
                        cbg_area, 
                        ndvi, 
                        tmax,
                        n, 
                        a_bb, 
                        fixmed, 
                        m_error),
                      nrow = 1)
        
        message(glue("Writing info for {ctf$species[i]} ind {ctf$ind_id[i]}, year {ctf$year[i]}, 
                   week {week}"))
        write.table(out, glue("{.outPF}/dbbmm_size.csv"), append = T, 
                    row.names = F, col.names = F, sep = ",")
        
      } #j
    }, error = function(e){cat(glue("ERROR: Size calulation failed for individual 
                                  {ctf$species[i]} {ctf$ind_id[i]}, year {ctf$year[i]}", 
                                  "\n"))})
    # if file hasn't been written
  }
}# if continue == T

#---- Finalize script ----#


message(glue('Script complete in {diffmin(t0)} minutes'))
