#!/usr/bin/env Rscript

# DESCRIPTION #
# 
# This script preps and cleans data for the COVID-19 Animal Movement Project
# See project documentation for details about anticipated directory structure.
# Major tasks of this script:
#   * Annotate dataset with the study periods:
#     * Pre-LD 2019
#     * LD 2019
#     * Post-LD 2019
#     * Pre-LD 2020
#     * LD 2020
#     * Post-LD 2020
#   * Remove events outside relevant study periods.
#   * Basic data cleaning:
#     * ID and remove outliers


'
Usage:
01-prep_and_clean.r [--db=<db>] 
01-prep_and_clean.r (-h | --help)

Control files:
  src/ctfs/dates.csv


Options:
-h --help     Show this screen.
-v --version     Show version.
-d --db=<db> Path to movement database. Defaults to <wd>/data/move.db
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/project/covid-19_movement'
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')

} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$db)

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
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

`%notin%` <- Negate(`%in%`)

#---- Load control files ----#
periods <- read_csv(file.path(.wd,'analysis/ctfs/dates.csv'),
                    col_types=list("date" = col_date(format = "%m/%d/%Y"))) 

#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

evt <- tbl(db,'event')
# dbBegin(db)
#---- Perform analysis ----#

#-- Make a filtered table by study period

# get list of inidviduals to remove from study base don bad coords
indtb <- tbl(db, "individual") %>%  collect()

rminds <- indtb %>% 
  filter(study_id == 1891172051 | study_id == 1891403240) %>% 
  pull(individual_id)

# extract only relevnt time periods
mod <- evt %>%
  filter((timestamp > !!periods$date[periods$cutpoint == "start_pre-ld_2019"] & 
            timestamp < !!periods$date[periods$cutpoint == "stop_2019"])
         | (timestamp > !!periods$date[periods$cutpoint == "start_pre-ld_2020"] & 
              timestamp < !!periods$date[periods$cutpoint == "stop_2020"])) %>%
  mutate(yr = strftime('%Y', timestamp),
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

#TODO: remove hardcode here
# filtering individual whose coords are utms not lat long
mod <- mod %>% 
  filter(individual_id %notin% rminds)

dbWriteTable(conn = db, name = "event_mod", value = mod, append = F, overwrite = T)

# # SQLite query to filter out events outside study periods, add a year variable,
# # and a study period variable.  Uses ctf to assign cutpoints
# q <- glue("
#   CREATE TABLE event_mod AS
#     SELECT *, strftime('%Y', `timestamp`) AS `yr`, 
#       CASE
#         WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_pre-ld_2019']}' 
#           AND `timestamp` < '{periods$date[periods$cutpoint == 'start_ld_2019']}') 
#             THEN ('pre-ld_2019')
#         WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_ld_2019']}' 
#           AND `timestamp` < '{periods$date[periods$cutpoint == 'start_post-ld_2019']}') 
#             THEN ('ld_2019')
#         WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_post-ld_2019']}' 
#           AND `timestamp` < '{periods$date[periods$cutpoint == 'stop_2019']}') 
#             THEN ('post-ld_2019')
#         WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_pre-ld_2020']}' 
#           AND `timestamp` < '{periods$date[periods$cutpoint == 'start_ld_2020']}') 
#             THEN ('pre-ld_2020')
#         WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_ld_2020']}' 
#           AND `timestamp` < '{periods$date[periods$cutpoint == 'start_post-ld_2020']}') 
#             THEN ('ld_2020')
#         WHEN (`timestamp` >= '{periods$date[periods$cutpoint == 'start_post-ld_2020']}' 
#           AND `timestamp` < '{periods$date[periods$cutpoint == 'stop_2020']}') 
#             THEN ('post-ld_2020')
#       END AS `trt`
#     FROM `event`
#       WHERE ((`timestamp` > '{periods$date[periods$cutpoint == 'start_pre-ld_2019']}' 
#         AND `timestamp` < '{periods$date[periods$cutpoint == 'stop_2019']}') 
#         OR (`timestamp` > '{periods$date[periods$cutpoint == 'start_pre-ld_2020']}' 
#         AND `timestamp` < '{periods$date[periods$cutpoint == 'stop_2020']}'))"
# )

# # Execute the query
# dbExecute(db, q)


#-- Clean Outliers

evt1 <- tbl(db, "event_mod")

cnt <- evt1 %>% 
  group_by(individual_id) %>% 
  summarise(n = n())

# get step lengths and turn angles across dataset
evt_sf <- evt1 %>% 
  # left_join(cnt) %>% 
  # filter(n < 3) %>% 
  as.data.frame() %>% 
  # group run calcs per individual
  group_by(individual_id) %>% 
  arrange(timestamp) %>% 
  mutate(rn = row_number(),
         # lead_d = geometry[row_number()+1],
         lag_lon = dplyr::lag(lon, 1),
         lag_lat = dplyr::lag(lat, 1),
         sl = distGeo(cbind(lon,lat), cbind(lag_lon, lag_lat)),
         bearing = bearing(cbind(lon,lat), cbind(lag_lon, lag_lat)),
         ta = 180-abs(180 - abs(bearing - dplyr::lag(bearing, 1)) %% 360))

cuts <- evt_sf %>% 
  as.data.frame() %>% 
  group_by(individual_id) %>% 
  summarize(
    qta = quantile(ta, probs = 0.95, na.rm = T, names = F),
    qsl = quantile(sl, probs = 0.95, na.rm = T, names = F)
  ) 

# filter outliers
out <- evt_sf %>% 
  # just join the cutpoints back to the dataset
  left_join(cuts) %>% 
  # conservative outlier thresh, must be past 95% quant for BOTH sl and TA
  filter(sl < qsl | ta < qta) %>% 
  ungroup()

dbWriteTable(conn = db, name = "event_clean", value = out, append = FALSE, overwrite = T)


#---- Finalize script ----#

dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
