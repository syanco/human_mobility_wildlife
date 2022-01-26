#!/usr/bin/env Rscript 
#
#CONDA: covid
#
# This script generates individual dynamic brownian bridge models and associated 
# UDs for migratory  periods.

# TODO:  The dBBMM paramaters (e.g., window size, margin, error, etc.) are 
# currently hardcoded.  Could be passed in as options to the script.
# TODO: verify the volume/probability problem for write out.
# 
# ==== Setup ====

'
Estimate Dynamic Brownian Bridge Movement Models (dBBMs) for pre-segemented 
migratory periods across multiple individuals

Usage:
make_dbbmm.r <db> <out> <nc>
make_dbbmm.r (-h | --help)

Parameters:
  db: path to movement databse. 
  out: path to output directory.
  nc: number of cores for parallel processing
  
Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/project/covid-19_movement'
  rd <- here::here
  
  .outPF <- file.path(.wd,'out')
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
  
  .nc <- 2
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  .nc <- ag$nc
  
}

#---- Initialize Environment ----#
# Start time
t0 <- Sys.time()

# Run startup
source(file.path(.wd,'analysis/src/startup.r'))

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
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))
indtb <- tbl(db,'individual') %>% 
  collect()

message("Disconnecting from databse...")
dbDisconnect(db)

ind <- indtb %>% 
  pull(individual_id)


yearvec <- c("2019", "2020")

registerDoMC(.nc)

# Toggle `%do%` to `%dopar%` for HPC, %do% for local
# foreach(j = 1:length(inds)) %do% {
foreach(j = 1:length(ind), .errorhandling = "pass", .inorder = F) %:%
  foreach(i = 1:2, .errorhandling = "pass", .inorder = F) %dopar% {
    message(glue("Starting individual {ind[j]}, year {yearvec[i]}..."))  
    
    #---- Initialize database ----#
    message(glue("Initializing database connection for individual {ind[j]}, year {yearvec[i]}..."))
    
    invisible(assert_that(file.exists(.dbPF)))
    db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
    invisible(assert_that(length(dbListTables(db))>0))
    
    
    #---- Perform analysis ----#
    message(glue("Gathering movement data for individual {ind[j]}, year {yearvec[i]}..."))
    
    evt0 <- tbl(db, "event_clean")
    
    scientificname <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(taxon_canonical_name)
    
    studyid <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(study_id)
    
    
    message("Filtering data and manipulating dates...")
    
    # prep data
    evt_mod <- evt0 %>% 
      # extract ind
      filter(individual_id == !!ind[j]) %>% 
      # extract year
      filter(yr == !!yearvec[i]) %>%
      collect() %>% 
      mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %T")) %>% 
      # sort by timestamp
      arrange(timestamp)
    
    dbDisconnect(db) 
    
    # TODO: this is an arbitrary minimum... check
    if(nrow(evt_mod) <= 25){
      message(glue("Insufficient records for individual {ind[j]}, year {yearvec[i]}, movin on..."))
      
      # Make entry in log file
      outlog <- matrix(c(scientificname, ind[j], studyid, yearvec[i], "dbbm", 
                         glue("dbbmm_{ind[j]}_{yearvec[i]}.rdata"),
                         0, as.character(Sys.Date())), 
                       nrow = 1)
      write.table(outlog, glue("{.outPF}/dbbmm_log.csv"), append = T, row.names = F, 
                  col.names = F, sep = ",")
      
    } else {
      
      #-- Fit dBBMMs
      message("Estimating dBBMMs...")
      tryCatch({
        # make minimal df for `move`    
        evt_tmp <- evt_mod %>% 
          select(lon, lat, timestamp, trt)
        
        evt_mv <- move(x=evt_tmp$lon, y=evt_tmp$lat, time=evt_tmp$timestamp, trt = evt_tmp$trt,
                       proj=CRS("+proj=longlat"))
        burstid <- factor(evt_tmp$trt[1:(n.locs(evt_mv)-1)])
        
        evt_burst <- burst(evt_mv, burstid)
        evt_mv_t <- spTransform(evt_burst, center = T)
        dbb_var <- brownian.motion.variance.dyn(object = evt_mv_t, 
                                                # TODO check on error modeling
                                                location.error = 10,
                                                #TODO: think about these...
                                                margin = 3, 
                                                window.size = 11)
        dbbm <- brownian.bridge.dyn(dbb_var, raster = 1000, 
                                    location.error = rep(10, n.locs(evt_mv_t)), 
                                    ext = 1,
                                    margin = 3, window.size = 11)
      }, error = function(e){cat("ERROR: unspecified error in fitting dBBMM", 
                                 "\n")})
      tryCatch({
        if(exists("dbbm")){
          # # Get UD volume
          # vol <- getVolumeUD(dbbm)
          # 
          # # Make 95% mask
          # mask95 <- vol
          # mask95[mask95>0.95] <- NA
          # mask95[mask95<0.95] <- 1
          # 
          # # Clip out 95% UD
          # ud95 <- dbbm*mask95
          
          # write the dbbmm objects to list
          tmp_out <- list("dBBMM Variance" = dbb_var,
                          "dBBMM Object" = dbbm,
                          # "Contours" =  raster2contour(dbbm),
                          # "UD Volume" = vol,
                          # "95% Mask" = mask95,
                          # "95% UD" = ud95,
                          "events" = evt_mod
          )
        } else {
          message(glue("No dbbmm object in memory, nothing written to tmp for individual {ind[j]}, {yearvec[i]}"))
        }
      }, error = function(e){cat("ERROR: couldnt write dBBMM objects to tmp", 
                                 "\n")})
      
      #-- Save individual output
      
      message(glue("Writing output for individual {ind[j]} to file..."))
      tryCatch({
        if(exists("tmp_out")){
          save(tmp_out,
               file = glue("{.outPF}/dbbmms/dbbmm_{ind[j]}_{yearvec[i]}.rdata")
          )
          
          # Make entry in log file
          outlog <- matrix(c(scientificname, ind[j], studyid, yearvec[i], "dbbm", 
                             glue("dbbmm_{ind[j]}_{yearvec[i]}.rdata"),
                             1, as.character(Sys.Date())), 
                           nrow = 1)
          write.table(outlog, glue("{.outPF}/dbbmm_log.csv"), append = T, row.names = F, 
                      col.names = F, sep = ",")
        }else {
          message(glue("No tmp list in memory, nothing written to file for {ind[j]}, {yearvec[i]}"))
        }
      }, error = function(e){cat("ERROR: couldnt save tmp_out to file", 
                                 "\n")})
      
    } # fi
  } #i (end loop through years) : #j (end loop through individuals)

#---- Finalize script ----#

message(glue('Script complete in {diffmin(t0)} minutes'))
