#!/usr/bin/env Rscript 
#
#CONDA: covid
#
# This script uses previously calculated dBBMMs to assess animal space use 
# during COVID-19 lockdowns

# TODO:  add migratory status filtering? Or maybe that's after...
# TODO: Update docopts
# TODO: Clean the script - there are unused object created, poorly commented, etc
# 
# ==== Setup ====

'
Calculate space use before and during COVID-19 lockdowns using previously estimated dBBMMs

Usage:
0X-calc_space_use.r <out> <db> <ctf> <nc>
0X-calc_space_use.r (-h | --help)

Parameters:
  --out: path to directory storing dBBMM outputs. 
  --db: path to dtabase
  --ctf: path to dbbmm log file
  --nc: number of cores for parallel processing
  
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
  .dbPF <- file.path(.wd,'processed_data/mosey_mod_20220303.db')
  .ctf <- file.path(.wd, "out/dbbmm_log.csv")
  
  .nc <- 4
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  .ctf <- makePath(ag$ctf)
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
    library(sf)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'), full.names=TRUE) %>%
  walk(source)


#---- Initialize database ----#
message("Initializing database connection and control files...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

evt_sg <- tbl(db, "event_sg") %>% 
  collect()

evt_cen <- tbl(db, "event_census") %>% 
  collect()

ind <- tbl(db,'individual') %>% 
  collect() %>% 
  pull(individual_id)

dbDisconnect(db)

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), file.path(.wd,'processed_data/mosey_mod_anno.db'), `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

evt_anno <- tbl(db, "event_clean") %>% 
  collect()

dbDisconnect(db)

# yearvec <- c("2019", "2020")
# trtvec <- c("pre-ld", "ld")



ctf <- read_csv(.ctf)
#TODO: rm below after a real run with log
ctf <- ctf[!duplicated(ctf),] %>% 
  filter(produced == 1)

# #- get the phylo control covar
# #TODO
# sp2019 <- raster("raw_data/geoserver_1643663843611/data.tif")
# sp2020 <- raster("raw_data/geoserver_1643663856740/data.tif")

#---- Perform analysis ----#
message("Gathering movement data...")



#+++++++++++++++++++++#
message(glue("Starting {.nc} cores"))
registerDoMC(.nc)

# Toggle `%do%` to `%dopar%` for HPC, %do% for local
foreach(i = 1:nrow(ctf), .errorhandling = "pass", .inorder = F) %dopar% {
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
      
      # # Get Phenology data
      # if(ctf$year[i] == "2020"){ #grab the correct phenology map
      #   # reproject the UD to match the spring data
      #   tmpr <- projectRaster(from = rb[[j]], to = sp2020)
      #   phen <- sum(values(tmpr*sp2020), na.rm = T)/ncell(tmpr[tmpr > 0])
      # } else {
      #   tmpr <- projectRaster(from = rb[[j]], to = sp2019)
      #   phen <- sum(values(tmpr*sp2019), na.rm = T)/ncell(tmpr[tmpr > 0])
      # } #else
      
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
        summarize(sg = mean(daily_count, na.rm = T))
      
      # Get pop density
      pop <- evt_cen %>%
        filter(event_id %in% evtids) %>% 
        summarize(pop = mean(total_population_2019, na.rm = T))
      
      # get Human encroachment
      # TODO: add code
      
      # Get NDVI
      ndvi <- evt_anno %>% 
        filter(event_id %in% evtids) %>% 
        summarize(sg = mean(ndvi, na.rm = T))
      
      # Get lst
      lst <- evt_anno %>% 
        filter(event_id %in% evtids) %>% 
        summarize(sg = mean(lst, na.rm = T))
      

      #unpack underlying data
      evt_tmp <- tmp_out$events 
      evt_mv <- move(x=evt_tmp$lon, y=evt_tmp$lat, time=evt_tmp$timestamp, trt = evt_tmp$trt,
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
      out <- matrix(c(ctf$species[i], ctf$ind_id[i], ctf$study_id[i], ctf$year[i], week, a, sg, pop, ndvi, lst, n, a_bb, fixmed, m_error),
                    nrow = 1)
      message(glue("Writing info for ind {ctf$ind_id[i]}, year {ctf$year[i]}, week {week}"))
      write.table(out, glue("{.outPF}/dbbmm_size.csv"), append = T, row.names = F, 
                  col.names = F, sep = ",")
      
    } #j
  }, error = function(e){cat(glue("ERROR: Size calulation failed for individual {ctf$ind_id[i]}, year {ctf$year[i]}", 
                                  "\n"))})
}

#---- Finalize script ----#


message(glue('Script complete in {diffmin(t0)} minutes'))
