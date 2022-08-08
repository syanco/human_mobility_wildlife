if(interactive()) {
  # rm(list=ls())
  library(here)
  
  .wd <- '/gpfs/ysm/project/jetz/ryo3/projects/covid'
  .test <- TRUE
  # rd <- here::here
  .datPF <- file.path(.wd,'analysis/')
  .outPF <- file.path(.wd,'analysis/')
  
  .dbPF <- '/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db'
  
} else {
  library(docopt)
  # library(rprojroot)
  
  .wd <- '/gpfs/loomis/project/jetz/sy522/covid-19_movement'
  # .script <-  thisfile()
  # rd <- is_rstudio_project$make_fix_file(.script)
  .datPF <- file.path(.wd,'out/')
  .outPF <- file.path(.wd,'out/')
  
  .dbPF <- '/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db'
  
}

source(file.path(.wd,'/src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(data.table)
    library(lubridate)
    library(sf)
    library(raster)
  }))

#---- Collect arguments ----#
args = commandArgs(trailingOnly = TRUE)

start_ix <- as.numeric(args[1])
end_ix <- as.numeric(args[2])
n <- as.numeric(args[3])

#---- Collect arguments ----#

files <- list.files(paste0(.outPF,'ssf-background-pts/individual-files'),full.names = TRUE)

evt <- data.table::rbindlist(lapply(files[start_ix:end_ix], data.table::fread)) %>%
  mutate("step_id" = c(1:nrow(.)),
         "date" = as.character(as_date(t2_))) %>%
  filter(abs(x2_) < 180) %>%
  filter(abs(y2_) < 90) %>%
  st_as_sf(coords = c("x2_","y2_"), crs="+proj=longlat +datum=WGS84", remove = FALSE) %>%
  st_make_valid()


message("reading in census block group geometries...")
cbg_sf <- st_read(paste0(.wd,"/raw_data/safegraph_open_census_data_2010_to_2019_geometry/cbg.geojson"))
cbg_area <- fread(paste0(.datPF, "event-annotations/cbg-area.csv"), colClasses = "character") %>%
  select(cbg_2010, cbg_area_m2)


message("running intersection with cbg geometries...")
evt_cbg <- st_intersection(evt,cbg_sf) %>%
  st_drop_geometry() %>%
  mutate(cbg_2010 = as.character(CensusBlockGroup)) %>%
  left_join(., cbg_area, by = "cbg_2010")


reformatted_files_daily <- list.files(paste0(.datPF,"safegraph/counties-dates-2-10-22-reformatted/daily-data"), full.names = TRUE)

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

head(evt_sg_ghm)

message("export individuals...")
ids <- unique(evt$individual_id)

for(i in 1:length(ids)){
  id <- ids[i]
  message(id)
  
  d <- evt_sg_ghm %>%
    filter(individual_id == id)
  
  fwrite(d, paste0(paste0(.outPF, "/ssf-background-points/annotated/individual-",id,".csv")))
}


message("done!")