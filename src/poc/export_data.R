# export data for Ben and Fede

db <- dbConnect(RSQLite::SQLite(), "projects/covid-19_movement/processed_data/mosey_mod.db", `synchronous` = NULL)

evnt_full <- tbl(db, "event_clean") %>% distinct() %>% collect() 
std_full  <- tbl(db, "study") %>% distinct() %>% collect()
ind_full  <- tbl(db, "individual") %>% distinct() %>% collect()

evnt_clean <- tbl(db, "event_clean") %>% distinct() %>% collect()

write_csv(evnt_full, file = "~/projects/covid-19_movement/out/data_export/event_full.csv")
write_csv(evnt_clean, file = "~/projects/covid-19_movement/out/data_export/event_clean.csv")
write_csv(std_full, file = "~/projects/covid-19_movement/out/data_export/study_full.csv")
write_csv(ind_full, file = "~/projects/covid-19_movement/out/data_export/individual_full.csv")