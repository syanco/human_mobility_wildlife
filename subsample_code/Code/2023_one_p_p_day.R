#!/usr/bin/env Rscript
#CONDA: covid
#---- Input Parameters ----#
library(here)
library(docopt)
library(rprojroot)

source(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/startup.r'))
source(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/funs/input_parse.r'))
#Source all files in the auto load funs directory
list.files(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/funs/auto'), full.names=T) %>% walk(source)
.dbPF <-'/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/Data/mosey_mod.db'
.wd <- '/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/'
.outPF <- paste0(.wd,'out/')

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
    require(MVNH)
    require(dplyr)
    require(doMPI)
  }))

#Source all files in the auto load funs directory
#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))
indtb <- tbl(db,'individual') %>%  # Load a tibble with all individual animals
  collect()

# Load the entire event table:
evt0 <- tbl(db, "event_clean")%>%
  collect()

evt0$timestamp_ymd = ymd_hms(evt0$timestamp)
evt0$jday = lubridate::yday(evt0$timestamp_ymd)
evt0$jday_year = paste0(evt0$jday, '-', evt0$yr)
evt0$individual_week = paste0(evt0$individual_id, '-',evt0$wk,'-', evt0$yr)

n_sample_per_week = evt0 %>%
  group_by(individual_id, yr, wk) %>%
  count() %>%
  ungroup()

summary(n_sample_per_week) # Meidan is 73


one_p_per_day = plyr::ddply(evt0, 'individual_id', function(x){
  plyr::ddply( x, 'jday_year', function(y){
    data.frame(
      who = unique(y$individual_id),
      jday_year = unique(y$jday_year),
      timestamp = y[1,]$timestamp_ymd,
      lon = y[1,]$lon,
      lat = y[1,]$lat,
      eventID = unique(y[1,]$event_id)
    )
  })
})



# save(one_p_per_day, file = '/gpfs/loomis/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/one_point_per_day.Rdata')
save(one_p_per_day, file = '/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/one_point_per_day.Rdata')



message("Disconnecting from databse...")
dbDisconnect(db)

ind <- indtb %>% 
  pull(individual_id)

ind = ind[ind  %in% unique(evt0$individual_id)]

# save(ind, file = '/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/indtb.Rdata')

event_sg <- read_csv('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/Data/event-annotation/event_sg.csv')
event_ghm <- read_csv('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/Data/event-annotation/event_ghm.csv')


one_p_per_day$event_id = one_p_per_day$eventID
one_p_per_day = one_p_per_day %>% select(-eventID)
class(one_p_per_day$event_id)
require(bit64)
event_sg$event_id = as.integer64(event_sg$event_id)
event_ghm$event_id = as.integer64(event_ghm$event_id)
# Subset one point per day:

one_p_per_day_sg = left_join(one_p_per_day, event_sg, by = 'event_id')
one_p_per_day_sg_ghm = left_join(one_p_per_day_sg, event_ghm, by = 'event_id')

one_p_per_day_sg_ghm = one_p_per_day_sg_ghm %>%
  mutate(timestamp = timestamp.x) %>% 
  dplyr::select(-timestamp.x, -timestamp.y, -date_hour, -StateFIPS, -CountyFIPS, -TractCode, -BlockGroup, -cbg_2010, -State, -County, -MTFCC, -cbg_area_m2)



  # cannot open compressed file '/gpfs/loomis/pi/jetz/de293/Anthropause/Anthropause_2023/out/one_point_per_day.Rdata', probable reason 'No such file or directory'



# save(one_p_per_day_sg_ghm, file = '/gpfs/loomis/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/one_point_per_day_human_anno.Rdata')
save(one_p_per_day_sg_ghm, file = '/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/one_point_per_day_human_anno.Rdata')

