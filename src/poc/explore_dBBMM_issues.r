library(tidyverse)
library(RSQLite)
library(DBI)
library(ggplot2)
library(lubridate)
library(sf)

size_sum <- read.csv("out/ind_area_04112022.csv")

.dbPF <- file.path('processed_data/mosey_mod.db')
db <- dbConnect(RSQLite::SQLite(), .dbPF)

evt <- tbl(db,'event')
evt_c <- tbl(db,'event_clean')
indtb <- tbl(db, "individual") %>%  collect()

ggplot(size_sum) +
  geom_boxplot(aes(x = sp, y = mean_area)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(size) +
  geom_boxplot(aes(x = species, y = area)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


targ_ind <- size_sum %>%
  filter(sp == "Scolopax_minor" | sp == "Rostrhamus_sociabilis") %>%
  # filter(sp == "Scolopax_minor") %>%
  pull(ind_f)

indtb %>%
  filter(individual_id %in% targ_ind) %>%
  pull(study_id)

x <- evt %>%
  filter(individual_id %in% targ_ind) %>%
  collect()

one <- x %>%
  filter(individual_id == 1891588271,
         timestamp > "2019-01-01")

one_clean <- evt_c %>%
  filter(individual_id == 1891588271) %>%
  collect()
####################################

evt_sf <- evt_tmp %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

ggplot() +
  geom_sf(data = evt_sf[evt_sf$wk == 1,], aes(color = as.factor(wk))) +
  coord_sf()
