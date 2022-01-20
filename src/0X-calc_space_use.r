#!/usr/bin/env Rscript 
#
#CONDA: covid
#
# This script uses previously calculated dBBMMs to assess animal space use 
# during COVID-19 lockdowns

# TODO:  add migratory status filtering
# TODO: Update docopts
# 
# ==== Setup ====

'
Calculate space use before and during COVID-19 lockdowns using previously estimated dBBMMs

Usage:
make_dbbmm.r <in> <out> <nc>
make_dbbmm.r (-h | --help)

Parameters:
  in: path to directory storing dBBMM outputs. 
  out: path to output directory.
  nc: number of cores for parallel processing
  
Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/covid-19_movement'
  rd <- here::here
  
  .outPF <- file.path(.wd,'out')
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
  
  .nc <- 2
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(rd('funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  .nc <- makePath(ag$nc)
  
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


#---- Perform analysis ----#
message("Gathering movement data...")

evt0 <- tbl(db, "event_clean")
indtb <- tbl(db,'individual')

yearvec <- c("2019", "2020")
trtvec <- c("pre-ld", "ld")

ind <- indtb %>% 
  collect() %>% 
  pull(individual_id)

registerDoMC(.nc)

# Toggle `%do%` to `%dopar%` for HPC, %do% for local
# foreach(j = 1:length(inds)) %do% {
foreach(j = 1:length(ind), i = 1:2, .errorhandling = "pass") %dopar% {
  message(glue("Starting individual {ind[j]}, year {yearvec[i]}..."))
  
  
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
    # sort by timestamp
    arrange(timestamp) %>% 
    collect()
  
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
    
    # make minimal df for `move`    
    evt_tmp <- evt_mod %>% 
      select(lon, lat, timestamp)
    
    evt_mv <- move(x=evt_tmp$lon, y=evt_tmp$lat, time=evt_tmp$timestamp, 
                   proj=CRS("+proj=longlat"))
    evt_mv_t <- spTransform(evt_mv, center = T)
    dbb_var <- brownian.motion.variance.dyn(object = evt_mv_t, 
                                            # TODO check on error modeling
                                            location.error = 1,
                                            #TODO: think about these...
                                            margin = 3, 
                                            window.size = 11)
    dbbm <- brownian.bridge.dyn(dbb_var_sp[[i]], raster = 1000, location.error = 1, 
                                margin = 3, window.size = 11)
    
    # Get UD volume
    vol <- getVolumeUD(dbbm)
    
    # Clip out 95% contour
    ud95 <- vol
    ud95[ud95>0.99] <- NA
    
    # write the dbbmm objects to list
    tmp_out <- list("dBBMM Variance" = dbb_var,
                    "dBBMM Object" = dbbm_sp,
                    "Contours" =  raster2contour(dbbm),
                    "UD Volume" = vol,
                    "95% Volume" = ud95,
                    "events" = evt_mod
    )
  } # fi
  
  #-- Save individual output
  
  # declare output destination
  .outTMP <- file.path(.outPF, scientificname, ind[j])
  
  message(glue("Writing output for individual {ind[j]} to file..."))
  
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
  
} #j (end loop through individuals)

#---- Finalize script ----#

message("Disconnecting from databse...")
dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
