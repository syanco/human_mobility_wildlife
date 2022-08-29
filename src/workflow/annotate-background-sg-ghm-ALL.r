if(interactive()) {
  # rm(list=ls())
  # library(here)
  
  .wd <- getwd()
  # .test <- TRUE
  # rd <- here::here
  .datPF <- file.path(.wd,'out/')
  .outPF <- file.path(.wd,'out/')
  
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
  .cores <- 1
  
} else {
  library(docopt)
  # library(rprojroot)
  
  .wd <- getwd()
  # .script <-  thisfile()
  # rd <- is_rstudio_project$make_fix_file(.script)
  .datPF <- file.path(.wd,'out/')
  .outPF <- file.path(.wd,'out/')
  
  .dbPF <- '/gpfs/loomis/pi/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db'
  
  #---- Collect arguments ----#
  args = commandArgs(trailingOnly = TRUE)
  
  .cores <- as.numeric(args[1])
}

source(file.path(.wd,'analysis/src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(data.table)
    library(lubridate)
    library(sf)
    library(raster)
    library(RSQLite)
    library(foreach)
    library(doMC)
  }))





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

# get ind count per species
sp_sum <- indtb %>%
  group_by(taxon_canonical_name) %>%
  summarize(nind = length(unique(individual_id))) 
# %>%
#   filter(nind > .minsp) #require a minimum of 10 individuals


# get focal species
sp <- sp_sum %>% 
  pull(taxon_canonical_name)



#---- Collect arguments ----#

registerDoMC(.cores)

foreach(i = 1:length(sp), .errorhandling = "pass", .inorder = F)%do%{
  
  # list of individuals within species
  indls <- indtb %>% 
    filter(taxon_canonical_name == sp[i]) %>% 
    pull(individual_id)
  
  foreach(j = 1:length(indls), .errorhandling = "pass", .inorder = F) %dopar% {
    f <- glue("out/ssf-background-pts/individual-files/{indls[j]}.csv")
    if(file.exists(f)){
      #check if we've already completed the annotation 
      #TURN THIS OFF TO OVERWRITE (OR TODO add an overwrite flag)
      if(file.exists( glue("{.outPF}ssf-background-pts/annotated/sindividual-{indls[j]}-sg-ghm.csv"))){
        evt <- read_csv(f)%>% 
          mutate("step_id" = c(1:nrow(.)),
                 "date" = as.character(as_date(t2_))) %>%
          filter(abs(x2_) < 180) %>%
          filter(abs(y2_) < 90) %>%
          st_as_sf(coords = c("x2_","y2_"), crs="+proj=longlat +datum=WGS84", remove = FALSE) %>%
          st_make_valid()
        
        message("reading in census block group geometries...")
        cbg_sf <- st_read(paste0(.wd,"/raw_data/safegraph_open_census_data_2010_to_2019_geometry/cbg.geojson"))
        cbg_area <- fread(paste0(.datPF, "event-annotation/cbg-area.csv"), colClasses = "character") %>%
          select(cbg_2010, cbg_area_m2)
        
        cbg_idx <- unique(unlist(st_intersects(evt, cbg_sf)))
        # reduc e the multipolygon object to only relevant features
        cbg_red <- cbg_sf[cbg_idx,]
        
        message("running intersection with cbg geometries...")
        # do the intersection on the reduced polygon
        evt_cbg <- st_intersection(evt, cbg_red) %>%
          st_drop_geometry() %>%
          mutate(cbg_2010 = as.character(CensusBlockGroup)) %>%
          left_join(., cbg_area, by = "cbg_2010")
        
        reformatted_files_daily <- list.files(paste0("processed_data/","safegraph/counties-dates-2-10-22-reformatted/daily-data"), full.names = TRUE)
        
        # combine all data
        message("reading in safegraph data...")
        daily_data <- data.table::rbindlist(lapply(reformatted_files_daily, data.table::fread, colClasses = "character"),use.names = TRUE) %>%
          select(cbg,date,count) %>%
          rename(safegraph_daily_count = count) %>%
          mutate(cbg_2010 = as.character(cbg),
                 date = as.character(date)) %>%
          select(-cbg)
        
        evt_sg <- left_join(evt_cbg,daily_data, by = c("cbg_2010", "date")) %>%
          select(-date)
        
        # read in global human modification layer
        message("reading in human modification...")
        ghm <- raster(paste0(.wd,"/raw_data/gHM/gHM.tif"))
        
        # transform to raster reference system
        message("transform event table...")
        evt_ghm <- st_transform(evt,st_crs(ghm))
        
        message("extract ghm...")
        evt_ghm$ghm <- raster::extract(ghm,evt_ghm)
        
        message("join annotations...")
        evt_ghm <- evt_ghm %>%
          st_drop_geometry()  %>%
          select(step_id, ghm)
        
        evt_sg_ghm <- left_join(evt_sg, evt_ghm, by = "step_id")
        
        # head(evt_sg_ghm)
        message(glue("exporting individual {j} of {length(indls)}..."))
        fwrite(evt_sg_ghm, glue("{.outPF}ssf-background-pts/annotated/sindividual-{indls[j]}-sg-ghm.csv"))
      }else{message(glue("{f} already annotated, skipping..."))}  
    } else {
      message(glue("{f} doesn't exist, moving on..."))
    } #else
  } #j
} #i


message("done!")