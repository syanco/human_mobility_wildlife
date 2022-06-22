library(sf)
library(RSQLite)
library(DBI)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)

.dbFP <- "processed_data/mosey_mod_20220303.db"

db <- dbConnect(RSQLite::SQLite(), .dbFP, 'synchronous' = NULL)

dat0 <- tbl(db, 'event_clean') %>%
  collect()

dbDisconnect(db)

dat <- dat0 %>% 
  mutate(
    ind_f = factor(individual_id),
    ts = date(timestamp),
    doy = yday(ts)) %>% 
  distinct(ind_f, yr, doy, .keep_all = T) 

dat_sf <- dat %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")

world <- ne_countries(returnclass = "sf")

NAm <- world %>% 
  filter(continent == "North America") %>% 
  st_transform("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")

NAbb <- st_bbox(dat_sf, n = 1)

dat_sf_filt <- dat_sf %>% 
  st_crop(NAbb)

ggplot() +
  geom_sf(data = NAm)+
  geom_sf(data = dat_sf_filt, color = "lightsalmon2", alpha = 0.2)+
  coord_sf(
    # xlim = c(-180, -20),
    # ylim = c(20, 80)
  ) +
  theme_minimal()+
  theme(axis.text = element_blank())
