#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# Trim database to only include studies from target period and area and
# correct species names in database. Generate new tables for the event, 
# study, and individual data that serve as input into the next workflow
# steps.

'
Trim data for the COVID-19 Animal Movement Project
- within the United States of America
- from specified periods within 2019-2020
- remove specific studies
- correct certain species names

See project documentation for details about anticipated directory structure.
Expects db to be of format mosey_db, with a table named "event"

Usage:
trim_data.r [--db=<db>] 
trim_data.r (-h | --help)

Control files:
  ctfs/dates.csv
  
Conda Environment: spatial

Options:
-h --help     Show this screen.
-v --version     Show version.
-d --db=<db> Path to movement database. Defaults to <wd>/data/move.db
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '/Users/juliet/Documents/OliverLab/covid_paper/repositories/human_mobility_wildlife'
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
  
} else {
  library(docopt)
  library(rprojroot)
  .wd <- getwd()
  source(file.path(.wd, 'src/funs/input_parse.r'))
  .dbPF <- '/tmp/mosey_mod.db'
  
}

#---- Initialize Environment ----#

# Start time
t0 <- Sys.time()

# Run startup
source(file.path(.wd,'src/startup.r'))

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(sf)
    library(geosphere)
    library(rnaturalearth)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)

`%notin%` <- Negate(`%in%`)

#---- Load control files ----#
periods <- read_csv(file.path(.wd,'ctfs/dates.csv'),
                    col_types=list("date" = col_date(format = "%m/%d/%Y")))

#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

evt <- tbl(db,'event') %>% collect()

#---- Perform analysis ----#

#-- Make a filtered table by study period

# Collect individual table
indtb <- tbl(db, "individual") %>%
         # correct test studies that have NA or just genus for taxon
         mutate(taxon_canonical_name = case_when(
        study_id == 2548691779 ~ "Odocoileus hemionus", # correct, keep
        study_id == 2575515057 ~ "Cervus elaphus", # correct, keep
        study_id == 1044238185 ~ "Alces alces", # correct, keep
        TRUE ~ taxon_canonical_name)) %>%
        collect()

# Collect study table
stdtb <- tbl(db, "study") %>% 
         collect() 


# extract only relevant time periods
mod <- evt %>%
  mutate(timestamp = as.POSIXct(timestamp, format="%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  filter((timestamp > !!periods$date[periods$cutpoint == "start_pre-ld_2019"] &
            timestamp < !!periods$date[periods$cutpoint == "stop_2019"])
          | (timestamp > !!periods$date[periods$cutpoint == "start_pre-ld_2020"] &
              timestamp < !!periods$date[periods$cutpoint == "stop_2020"])) %>%
  mutate(yr = strftime(timestamp, '%Y'),
         trt = case_when(
           timestamp >= !!periods$date[periods$cutpoint == "start_pre-ld_2019"] &
             timestamp < !!periods$date[periods$cutpoint == "start_ld_2019"] ~ "pre-ld_2019",
           timestamp >= !!periods$date[periods$cutpoint == "start_ld_2019"] &
             timestamp < !!periods$date[periods$cutpoint == "start_post-ld_2019"] ~ "ld_2019",
           timestamp >= !!periods$date[periods$cutpoint == "start_post-ld_2019"] &
             timestamp < !!periods$date[periods$cutpoint == "stop_2019"] ~ "post-ld_2019",
           timestamp >= !!periods$date[periods$cutpoint == "start_pre-ld_2020"] &
             timestamp < !!periods$date[periods$cutpoint == "start_ld_2020"] ~ "pre-ld_2020",
           timestamp >= !!periods$date[periods$cutpoint == "start_ld_2020"] &
             timestamp < !!periods$date[periods$cutpoint == "start_post-ld_2020"] ~ "ld_2020",
           timestamp >= !!periods$date[periods$cutpoint == "start_post-ld_2020"] &
             timestamp < !!periods$date[periods$cutpoint == "stop_2020"] ~ "post-ld_2020")
  ) %>%
  collect()

mod2 <- mod %>%
  left_join(indtb %>% select(-study_id), # drop study_id col in ind data before joining
            by = "individual_id") %>%
  mutate(ind_f = as.factor(individual_id),
         species = taxon_canonical_name)%>%  # create factor version of ind for REs)
  mutate(species = case_when( # correct species names
    study_id == 1442516400 ~ "Anser caerulescens", # correct, keep
    study_id == 1631574074 ~ "Ursus americanus", # correct, keep
    study_id == 1418296656 ~ "Numenius americanus", # correct, keep
    study_id == 474651680  ~ "Odocoileus virginianus", # correct, keep
    study_id == 1044238185 ~ "Alces alces", # correct, keep
    study_id == 2548691779 ~ "Odocoileus hemionus", # correct, keep
    study_id == 2575515057 ~ "Cervus elaphus", # correct, keep
    TRUE ~ species
  ))%>%
  mutate(species = case_when(
    species == "Chen caerulescens" ~ "Anser caerulescens", # correct, keep
    species == "Chen rossii" ~ "Anser rossii", # correct, keep
    TRUE ~ species
  ))

# Filter to US only #
mod2_sf <- st_as_sf(mod2, coords = c("lon", "lat"), crs = 4326)

us <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name == "United States of America")

# Find cells that overlap land in US
keep_fixes <- apply(st_intersects(mod2_sf, us, sparse = F),
                    1,
                    function(x){
                      ifelse(sum(x[1])>0, T, F)
                    }
)

# Filter grid to only cells on land
mod3 <- mod2[keep_fixes,]

# convert the timestamp class of POSIXct back to character bc SQLite cannot handle class POSIXct
mod3 <- mod3 %>%
  mutate(timestamp = format(timestamp, "%Y-%m-%d %H:%M:%S"))

# write results back to db
dbWriteTable(conn = db, name = "event_trim", value = mod3, append = F, overwrite = T)

# Correct the individual table
ind_out <- indtb %>%
  semi_join(mod3, by = "individual_id")   # only retain inds contained in cleaned event table

# write table back to db
dbWriteTable(conn = db, name = "individual_trim", value = ind_out, append = FALSE, overwrite = T)


# Correct the study table
std_out <- stdtb %>%
  semi_join(ind_out, by = "study_id") # only retain studies contained in cleaned individual table

# write table back to db
dbWriteTable(conn = db, name = "study_trim", value = std_out, append = FALSE, overwrite = T)


#---- Finalize script ----#

# disconnect from db
dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
