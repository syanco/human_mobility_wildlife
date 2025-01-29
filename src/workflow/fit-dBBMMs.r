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
  
  .wd <- getwd()
  rd <- here::here
  
  .outPF <- file.path(.wd,'out')
  .dbPF <- file.path('/Users/juliet/Documents/OliverLab/covid_paper/db_explore/db_after_hmw_wf_part1/mosey_mod_clean-movement_complete.db')
  
  .nc <- 2
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$db)
  .outPF <- makePath(ag$out)
  .nc <- ag$nc
  
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
    library(glue)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))
indtb <- tbl(db,'individual_clean') %>% 
  collect() 

indtb <- indtb[!duplicated(indtb),]

message("Disconnecting from database...")
dbDisconnect(db)

ind <- indtb %>% 
  pull(individual_id)

#remove duplicated created due to merging of dbs
ind <- ind[!duplicated(ind)] 

yearvec <- c("2019", "2020")

# Read in output log
log <- read_csv(glue("{.outPF}/dbbmm_log.csv"))

registerDoMC(.nc)

# Toggle `%do%` to `%dopar%` for HPC, %do% for local
# foreach(j = 1:length(inds)) %do% {
foreach(j = 1:length(ind), .errorhandling = "pass", .inorder = F) %:%
  foreach(i = 1:2, .errorhandling = "pass", .inorder = F) %dopar% {
    # check whether individual has been previously considered
    if(log %>% 
       filter(ind_id == ind[j] & year == yearvec[i]) %>% 
       nrow() == 0) {
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
        distinct() %>% 
        # sort by timestamp
        arrange(timestamp)
      
      dbDisconnect(db)
      
      # check that the filtered dataframe has any data
      if (nrow(evt_mod) == 0) {
        
        # if no rows, record this combination as failed in the appropriate df
        # in order to track how many errors are from lack of data
        # versus failed dBBMM generation
        no_ind_yr_pairs <- data.frame("study_id" = studyid,
                                      "individual_id" = ind[j],
                                      "year" = yearvec[i])
        
        write.table(no_ind_yr_pairs, 
                    glue("{.outPF}/no_ind_yr_pairs.csv"), 
                    append = T, 
                    row.names = F, 
                    col.names = F, 
                    sep = ",")

        message(glue("No paired data available for ind {ind[j]}, yr {yearvec[i]}; no dBBMM available."))
        
        # move onto the next iteration of ind and year, no use in trying to fit the dMMBB
        return(NULL)
        
      } else {
        
        # if event data is large...
        # TODO:  I think this upper size check is unecessary - maybe remove someday, but for now I just set the threshold very high
        if(nrow(evt_mod) > 100000){
          
          # ...make an entry in the big mem log file - save dbbmm for later
          outlog <- data.frame("species" = scientificname, 
                               "ind_id" = ind[j], 
                               "study_id" = studyid, 
                               "year" = yearvec[i], 
                               "out_type" = "dbbm", 
                               "filename" = glue("dbbmm_{ind[j]}_{yearvec[i]}.rdata"), 
                               "produced" = 0, 
                               "out_date" = as.character(Sys.Date()))
          write.table(outlog, glue("{.outPF}/dbbmm_bigmem_log.csv"), append = T, row.names = F, 
                      col.names = F, sep = ",")
        } else {
          # ...else proceed with fitting dbbmm
          
          
          #-- Fit dBBMMs
          message("Estimating dBBMMs...")
          tryCatch({
            # make minimal df for `move`    
            evt_tmp <- evt_mod %>% 
              select(lon, lat, timestamp, wk)
            
            evt_mv <- move(x=evt_tmp$lon, y=evt_tmp$lat, time=evt_tmp$timestamp,
                           proj=CRS("+proj=longlat"))
            burstid <- factor(evt_tmp$wk[1:(n.locs(evt_mv)-1)])
            #id "intended fix rate"
            fixmed <- median(timeLag(x=evt_mv, units="mins"))
            evt_burst <- burst(evt_mv, burstid)

            # transform for spTransform's default CRS, equidistant
            # evt_mv_t <- spTransform(evt_burst, center = T)
            # transform to an equal area CRS (Mollweide), as is advised by documentation for dbbmm documentation
            evt_mv_t <- spTransform(evt_burst, CRSobj="+proj=moll +ellps=WGS84")

            # remove variance of the segments corresponding to large time gaps
            dbb_var <- brownian.motion.variance.dyn(object = evt_mv_t, 
                                                    location.error = if(any(is.na(evt_mod$horizontal_accuracy))){
                                                      rep(5, n.locs(evt_mv_t))}else{
                                                        evt_mod$horizontal_accuracy},
                                                    margin = 11, 
                                                    window.size = 31)
            # remove any segments with gaps > 3x the intended fix rate
            dbb_var@interest[timeLag(evt_mv_t,"mins")>(fixmed*3)] <- FALSE
            
            dbbm <- brownian.bridge.dyn(dbb_var, 
                                        location.error = if(any(is.na(evt_mod$horizontal_accuracy))){
                                          rep(5, n.locs(evt_mv_t))}else{
                                            evt_mod$horizontal_accuracy} ,
                                        time.step = (fixmed/15),
                                        dimSize = 1000,
                                        ext = 10,
                                        margin = 11, window.size = 31)
          }, error = function(e){cat(glue("ERROR: unspecified error in fitting dBBMM for ind {ind[j]}, yr {yearvec[i]}", 
                                          "\n"))})
          tryCatch({
            
            if(exists("dbbm")){
              # write the dbbmm objects to list
              tmp_out <- list("dBBMM Variance" = dbb_var,
                              "dBBMM Object" = dbbm,
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
              outlog <- data.frame("species" = scientificname, 
                                   "ind_id" = ind[j], 
                                   "study_id" = studyid, 
                                   "year" = yearvec[i], 
                                   "out_type" = "dbbm", 
                                   "filename" = glue("dbbmm_{ind[j]}_{yearvec[i]}.rdata"), 
                                   "produced" = 1, 
                                   "out_date" = as.character(Sys.Date()))
              write.table(outlog, glue("{.outPF}/dbbmm_log.csv"), append = T, row.names = F, 
                          col.names = F, sep = ",")
            } else {
              message(glue("No tmp list in memory, nothing written to file for {ind[j]}, {yearvec[i]}"))
            }
          }, error = function(e){cat("ERROR: couldnt save tmp_out to file", 
                                     "\n")})
        } # fi data not too large
      } # fi end the check whether the filtered dataframe has any data
    } # fi end the check whether individual has been previously considered
  } #i (end loop through years) : #j (end loop through individuals)

#---- Finalize script ----#

message(glue('Script complete in {diffmin(t0)} minutes'))
