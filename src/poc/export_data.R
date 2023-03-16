# export data for Ben and Fede
library(RSQLite)
library(DBI)
library(tidyverse)


db <- dbConnect(RSQLite::SQLite(), "~/projects/covid-19_movement/processed_data/mosey_mod.db", `synchronous` = NULL)

evnt_full <- tbl(db, "event") %>% distinct() %>% collect() 
write_csv(evnt_full, file = "~/projects/covid-19_movement/out/data_export/event_full.csv")

std_full  <- tbl(db, "study") %>% distinct() %>% collect()
write_csv(std_full, file = "~/projects/covid-19_movement/out/data_export/study_full.csv")

ind_full  <- tbl(db, "individual") %>% distinct() %>% collect()
write_csv(ind_full, file = "~/projects/covid-19_movement/out/data_export/individual_full.csv")

evnt_clean <- tbl(db, "event_clean") %>% distinct() %>% collect()
write_csv(evnt_clean, file = "~/projects/covid-19_movement/out/data_export/event_clean.csv")
