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
make_dbbmm.r <in> <out> <nc> <ctf>
make_dbbmm.r (-h | --help)

Parameters:
  in: path to directory storing dBBMM outputs. 
  out: path to output directory.
  nc: number of cores for parallel processing
  ctf: path to dbbmm log file
  
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
  .ctf <- file.path(.wd, "out/dbbmm_log.csv")
  
  .nc <- 2
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(rd('funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  ctf <- makePath(ag$ctf)
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
message("Initializing database connection and control files...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

ctf <- read_csv(.ctf)
#TODO: rm below after a real run with log
ctf <- ctf[!duplicated(ctf),] %>% 
  filter(produced == 1)


#---- Perform analysis ----#
message("Gathering movement data...")

evt0 <- tbl(db, "event_clean")
indtb <- tbl(db,'individual')

yearvec <- c("2019", "2020")
trtvec <- c("pre-ld", "ld")

ind <- indtb %>% 
  collect() %>% 
  pull(individual_id)

#+++++++++++++++++++++#
for(i in 1:nrow(ctf)){
  message(glue("Starting ind {ctf$ind_id[i]}, year {ctf$year[i]}"))
  tryCatch({
  load(glue("{.outPF}/dbbmms/dbbmm_{ctf$ind_id[i]}_{ctf$year[i]}.rdata"))
  
  r <- tmp_out$`dBBMM Object`
  # r
  # plot(sqrt(r))
  rb <- UDStack(r)
  UDr <- getVolumeUD(rb)
  # plot(UDr)
  for(j in 1:nlayers(UDr)){
    a <- ncell(UDr[[j]])/1000
    trt <- names(UDr[[j]])
    out <- matrix(c(ctf$species[i], ctf$ind_id[i], ctf$study_id[i], ctf$year[i], trt = trt, a),
                  nrow = 1)
    message(glue("Writing info for ind {ctf$ind_id[i]}, year {ctf$year[i]}, period {trt}"))
    write.table(out, glue("{.outPF}/dbbmm_size.csv"), append = T, row.names = F, 
                col.names = F, sep = ",")
    
  }
  }, error = function(e){cat(glue("ERROR: Size calulation failed for individual {ctf$ind_id[i]}, year {ctf$year[i]}", 
                            "\n"))})
}

#---- Finalize script ----#

message("Disconnecting from databse...")
dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
