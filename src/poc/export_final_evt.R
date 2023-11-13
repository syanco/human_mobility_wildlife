

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- getwd()
  rd <- here::here
  
  .outPF <- file.path(.wd,'out')
  .dbPF <- file.path(.wd,'processed_data/mosey_mod_2023.db')
  
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
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))
indtb <- tbl(db,'individual_clean') %>% 
  collect() 

evttb <- tbl(db, "event_clean") %>% 
  mutate(taxon_canonical_name = case_when( # correct species names
    study_id == 1442516400 ~ "Anser caerulescens",
    study_id == 1233029719 ~ "Odocoileus virginianus",
    study_id == 1631574074 ~ "Ursus americanus",
    study_id == 1418296656 ~ "Numenius americanus",
    study_id == 474651680  ~ "Odocoileus virginianus",
    study_id == 1044238185 ~ "Alces alces",
    TRUE ~ taxon_canonical_name
  ))%>% 
  mutate(taxon_canonical_name = case_when(
    taxon_canonical_name == "Chen caerulescens" ~ "Anser caerulescens",
    TRUE ~ taxon_canonical_name
  )) %>% 
  filter(taxon_canonical_name != "Rallus longirostris" & 
           taxon_canonical_name != "Numenius americanus" &
           taxon_canonical_name != "Buteo regalis" &
           taxon_canonical_name != "Cathartes aura" &
           taxon_canonical_name != "Canis lupus"
         ) %>% 
  collect()

length(unique(evttb$taxon_canonical_name))
sp_sum <- evttb %>%
  group_by(taxon_canonical_name) %>%
  summarize(nind = length(unique(individual_id)))
dbDisconnect(db)

write_csv(evttb, "out/final_evt_for_map.csv")
