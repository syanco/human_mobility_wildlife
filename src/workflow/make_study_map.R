#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  # .wd <- '~/projects/covid-19_movement'
  .wd <- '~/Documents/covid-19_movement/'
  
  .dbPF <- file.path(.wd,'processed_data/mosey_mod_2023.db')
  .wkmin <- 30
  .minsp <- 5
  
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$db)
  .wkmin <- ag$wkmin
  .minsp <- ag$minsp
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
    library(sf)
    library(geosphere)
    library(lubridate)
    library(rnaturalearth)
    library(ggsflabel)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

`%notin%` <- Negate(`%in%`)


#---- Initialize database ----#
message("Connecting to database...")
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

#---- Perform analysis ----#

# dbListTables(db)

#-- Load Cleaned Data
evt0 <- tbl(db, 'event_clean') %>% collect()
study_vec <- unique(evt_cln$study_id)

std0 <- tbl(db, 'study_clean') %>% 
  collect() %>% 
  filter(study_id %in% study_vec) %>% 
  mutate(species = case_when(study_id == 619097045 ~ "USGS Ducks",
                             TRUE ~ taxon_ids))

unique(std0$species)

std_sf <- st_as_sf(std0, coords = c("main_location_long", "main_location_lat"),
                   crs = 4326)

us <- ne_countries(country = "United States of America", returnclass = "sf")

ggplot()+
  geom_sf(data = us)+
  geom_sf(data = std_sf, inherit.aes = F)+
  # geom_sf_label_repel(data = std_sf, aes(label = species), force_pull = 0, force = 10) +
  coord_sf() + 
  # theme(legend.position = "bottom") +
  facet_wrap(~species)
  NULL

save(std_sf, file = "out/study_locs.rdata")
