# This script contains code to check on sample sizes prior to any modeling steps 
# (to figure out how and where species are getting filtered out of the anlysis)

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- getwd()
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
`%notin%` <- Negate(`%in%`)

# Run startup
source(file.path(.wd,'analysis/src/startup.r'))

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
    library(raster)
    # library(move)
    # library(doMC)
    library(foreach)
    library(mapview)
    library(sf)
    library(tidyverse)
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

indtb <- indtb[!duplicated(indtb),]

evttb <- tbl(db, 'event_clean') %>% 
  collect() %>% 
  left_join(indtb, by = c("individual_id", "study_id"))

message("Disconnecting from databse...")
dbDisconnect(db)


evt1 <- evttb %>% 
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>% 
  mutate(ind_f = as.factor(individual_id),
         species = taxon_canonical_name)%>%  # create factor version of ind for REs)
  mutate(species = case_when( # correct species names
    study_id == 1442516400 ~ "Anser caerulescens",
    study_id == 1233029719 ~ "Odocoileus virginianus",
    study_id == 1631574074 ~ "Ursus americanus",
    study_id == 1418296656 ~ "Numenius americanus",
    study_id == 474651680  ~ "Odocoileus virginianus",
    study_id == 1044238185 ~ "Alces alces",
    study_id == 2548691779 ~ "Odocoileus hemionus", # !!! CORRECT IN WORKFLOW !!! #
    study_id == 2575515057 ~ "Cervus elaphus",  # !!! CORRECT IN WORKFLOW !!! #
    TRUE ~ species
  ))%>% 
  mutate(species = case_when(
    species == "Chen caerulescens" ~ "Anser caerulescens",
    TRUE ~ species
  )) %>% 
  distinct() 


# How many species in data?
evt1 %>% 
  pull(species) %>% 
  unique() %>% 
  length()

# How many have any complete cases?
evt1 %>% 
  select(species, tmax, ndvi) %>% 
  filter(complete.cases(.)) %>% 
  pull(species) %>% 
  unique() %>% 
  length()


# List with at least some complete cases
(complete_spp <- evt1 %>% 
    select(species, tmax, ndvi) %>% 
    filter(complete.cases(.)) %>% 
    pull(species) %>% 
    unique())

# List with all spp.
(all_spp <- evt1 %>% 
    pull(species) %>% 
    unique() %>% 
    length())

# For which species did annotation totally fail
evt1 %>% 
  filter(species %notin% complete_spp) %>% 
  pull(species) %>% 
  unique()

# Individual ids of those that failed
(ind_fails <- evt1 %>% 
    filter(species %notin% complete_spp) %>% 
    pull(individual_id) %>% 
    unique() %>% 
    as.integer())




# Full data for non annotated
failed_dat <- evt1 %>% 
  filter(species %notin% complete_spp)

summary(failed_dat$ndvi)
summary(failed_dat$tmax)

# number of indibviduals for species if the annos hadn't failed...
# by spp
failed_dat %>% 
  group_by(species) %>% 
  summarize(n = n_distinct(ind_f))

# total
failed_dat %>% 
  summarize(n = n_distinct(ind_f))

# What studies
(std_fail <- unique(failed_dat$study_id))

# Make new study.csv file with only the faled individuals
std_ctf <- read.csv("analysis/ctfs/study.csv")
ctf_out <- std_ctf %>% 
  mutate(run = case_when(study_id == std_fail ~ 1,
                         TRUE ~ 0))
write_csv(ctf_out, file = "analysis/ctfs/study.csv")

# Convert to sf format
failed_sf <- st_as_sf(failed_dat, coords = c("lon", "lat"), crs = 4326)

# Quick map of failed
mapview(failed_sf)

# sample sizes of complete cases:
print(evt1 %>% 
        select(species, tmax, ndvi) %>% 
        filter(complete.cases(.)) %>% 
        group_by(species) %>% 
        summarise(n = n()), n= 100)


# print n unique inds:
print(evt1 %>% 
        select(species, tmax, ndvi, ind_f) %>% 
        filter(complete.cases(.)) %>% 
        group_by(species) %>% 
        summarise(n = n_distinct(ind_f)), n= 100)
