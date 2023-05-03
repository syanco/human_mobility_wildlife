#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- getwd()
  rd <- here::here
  
  .outP <- file.path(.wd,'out')
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
  
  .nc <- 2
  
} else {
  # library(docopt)
  # library(rprojroot)
  
  # ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .outP <- file.path(.wd, 'out')
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
  # .nc <- ag$nc
  
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
    library(glue)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

# load db tables
inddb <- tbl(db,'individual') 
evttb <- tbl(db, 'event_clean')
# indtb <- tbl(db,'individual') %>% 
#   collect() 

# load ghm, cbg, and sg
ghm <- read_csv("out/event-annotation/event_ghm.csv")
sg <- read_csv("out/event-annotation/event_sg.csv")


#----  Make complete annotated dataset ----#

evt_full <- evttb %>% 
  # head(n = 25) %>% 
  left_join(inddb) %>% 
  collect() %>% 
  mutate(event_id = as.numeric(event_id),
         timestamp = ymd_hms(timestamp)) %>% 
  mutate(ID = as.factor(individual_id)) %>% 
  left_join(ghm, by = c("event_id")) %>% 
  left_join(sg, by = c("event_id", "timestamp"))

saveRDS(evt_full, file = glue("{.outP}/anno_full_{Sys.Date()}.rds"))


#---- Make daily data ----#

evt_day <- evt_full %>% 
  mutate(doy = yday(timestamp)) %>% 
  group_by(ID, doy) %>% 
  summarize(species = taxon_canonical_name[1],
            ghm = mean(ghm, na.rm = T),
            sg = mean(safegraph_daily_count, na.rm = T),
            cbg = mean(cbg_area_m2, na.rm = T),
            sg_norm = sg/cbg,
            tmax = mean(tmax, na.rm = T),
            ndvi = mean(ndvi, na.rm = T),
            elev = mean(elev, na.rm = T),
            ta = mean(ta, na.rm = T),
            bearing = mean(bearing, na.rm = T),
            sl = mean(sl, na.rm = T),
            v = mean(v, na.rm = T),
            lat = mean(lat, na.rm = T),
            lon = mean(lon, na.rm = T)) %>% 
  ungroup()

saveRDS(evt_day, file = glue("{.outP}/anno_day_{Sys.Date()}.rds"))


#---- Make weekly data ----#

evt_wk <- evt_full %>% 
  group_by(ID, wk) %>% 
  summarize(species = taxon_canonical_name[1],
            ghm = mean(ghm, na.rm = T),
            sg = mean(safegraph_daily_count, na.rm = T),
            cbg = mean(cbg_area_m2, na.rm = T),
            sg_norm = sg/cbg,
            tmax = mean(tmax, na.rm = T),
            ndvi = mean(ndvi, na.rm = T),
            elev = mean(elev, na.rm = T),
            ta = mean(ta, na.rm = T),
            bearing = mean(bearing, na.rm = T),
            sl = mean(sl, na.rm = T),
            v = mean(v, na.rm = T),
            lat = mean(lat, na.rm = T),
            lon = mean(lon, na.rm = T)) %>% 
  ungroup()

saveRDS(evt_wk, file = glue("{.outP}/anno_week_{Sys.Date()}.rds"))


#---- Make weekly data ----#

evt_mo <- evt_full %>% 
  mutate(month = month(timestamp)) %>% 
  group_by(ID, month) %>% 
  summarize(species = taxon_canonical_name[1],
            ghm = mean(ghm, na.rm = T),
            sg = mean(safegraph_daily_count, na.rm = T),
            cbg = mean(cbg_area_m2, na.rm = T),
            sg_norm = sg/cbg,
            tmax = mean(tmax, na.rm = T),
            ndvi = mean(ndvi, na.rm = T),
            elev = mean(elev, na.rm = T),
            ta = mean(ta, na.rm = T),
            bearing = mean(bearing, na.rm = T),
            sl = mean(sl, na.rm = T),
            v = mean(v, na.rm = T),
            lat = mean(lat, na.rm = T),
            lon = mean(lon, na.rm = T)) %>% 
  ungroup()

saveRDS(evt_mo, file = glue("{.outP}/anno_month_{Sys.Date()}.rds"))