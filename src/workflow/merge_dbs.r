#!/usr/bin/env Rscript 
#
#CONDA: covid
#
# This scriptMerges the new and replacement studies in mosey_swap_mod.db into the 
# main analysis db (mosey_mod.db)

'
Merge mosey_swap_mod.db with mosey_mod.db

Usage:
merge_dbs.r 
merge_dbs.r (-h | --help)

Parameters:
  db_swap: path to the swap databse. 
  db_main: path to the main database.

Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/covid-19_movement'
  rd <- here::here
  
  # .dbmain <- '~/project/covid-19_movement/processed_data/mosey_mod.db'
  .dbmain <- '~/projects/covid-19_movement/processed_data/mosey_mod.db'
  .dbswap <- '~/projects/covid-19_movement/processed_data/mosey_swap.db'
  # .dbloomis <- '~/projects/covid-19_movement/processed_data/mosey_loomis.db'
  
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  # .dbmain <- '~/project/covid-19_movement/processed_data/mosey_mod.db'
  .dbmain <- '~/project/covid-19_movement/processed_data/mosey_mod.db'
  .dbswap <- '~/project/covid-19_movement/processed_data/mosey_swap.db'

  
}

#---- Initialize Environment ----#
# Start time
t0 <- Sys.time()

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
    library(tidyverse)
  }))

#---- Load Control Files ----#

# control file with new studies to add to the db
ctfnew <- read.csv('~/projects/covid-19_movement/ctfs/new_studies.csv', 
                   colClasses = c(Study_id = "character"))

# control file with studies to swap
ctfswap <- read.csv('~/projects/covid-19_movement/ctfs/swaps.csv', 
                    colClasses = c(Old_study_id = "character",
                                   New_study_id = "character"))
  
swapout <- ctfswap$Old_study_id  # vec of studies to swap out (not including the new ones)
# swapout_names <- ctfswap$Old_study_name
swapin <- ctfswap$New_study_id

newstudies <- ctfnew$Study_id

# Run startup
source(file.path(.wd,'analysis/src/startup.r'))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

# Send merge log to file
sink("out/merge_log.txt", split = T, append = F)

#---- Initialize database ----#
cat("Initializing database connection... \n")

invisible(assert_that(file.exists(.dbmain)))
dbmain <- dbConnect(RSQLite::SQLite(), .dbmain, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(dbmain))>0))

invisible(assert_that(file.exists(.dbswap)))
dbswap <- dbConnect(RSQLite::SQLite(), .dbswap, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(dbswap))>0))

# invisible(assert_that(file.exists(.dbswap)))
# dbloomis <- dbConnect(RSQLite::SQLite(), "processed_data/backups/mosey_loomis.db", `synchronous` = NULL)
# invisible(assert_that(length(dbListTables(dbswap))>0))

#---- Merge Study Tables ----3

cat("\n Merging study tables... \n")

std_swap <- tbl(dbswap, 'study') %>% 
  mutate(study_id = as.character(study_id)) %>% 
  collect()

std_main <- tbl(dbmain,'study') %>% # load main db event table
  mutate(study_id = as.character(study_id)) %>% 
  collect() 
# std_loomis <- tbl(dbloomis, 'study') %>% 
#   mutate(study_id = as.character(study_id)) %>% 
#   collect()

#check that colnames match
cat(glue("Do all column headings match? {all(colnames(std_swap) == colnames(std_main))}"), " \n")

# Pre swap counts and checks
cat(glue("Number of studies in db to be deleted: {sum(std_main$study_id %in% swapout)} \n"), " \n")
cat(glue("Number of studies to be replaced: {sum(std_swap$study_id %in% swapin)} \n"), " \n")
cat(glue("Did the swap change the count of studies? {sum(std_main$study_id %in% swapout) != sum(std_swap$study_id %in% swapin)} \n"), " \n")
cat(glue("No replacement studies currently in individual table? {(sum(std_main$study_id %in% swapin))==0} \n"), " \n")

cat(glue("Number of new studies to add to db: {sum(std_swap$study_id %in% newstudies)} \n"), " \n")
cat(glue("No new studies currently in db? {(sum(std_main$study_id %in% newstudies))==0} \n"), " \n")

std_final_pred <- nrow(std_main) - sum(std_main$study_id %in% swapout) + sum(std_swap$study_id %in% swapin) + sum(std_swap$study_id %in% newstudies)
cat(glue("Final number of studies should be: {std_final_pred} \n"), " \n")

std_out <- std_main %>%  
  filter(study_id %notin% swapout) %>% # remove studies to be swapped
  bind_rows(std_swap)   # bind new event table in

cat(glue("\n Final table length as expected? {nrow(std_out)==std_final_pred} \n"), " \n")

# write back to db
cat("Writing study table back to db... \n ")
dbWriteTable(dbmain, 'study', std_out, overwrite = T)

#---- Merge Individual Tables ----#

cat("\n Merging individual tables... \n")

ind_swap <- tbl(dbswap, 'individual') %>% # load the swap ind table
  mutate(study_id = as.character(study_id)) %>% 
  collect()

ind_main <- tbl(dbmain,'individual') %>% # load main db ind table
  mutate(study_id = as.character(study_id)) %>% 
  collect() 

# ind_loomis <- tbl(dbloomis,'individual') %>% # load main db ind table
#   mutate(study_id = as.character(study_id)) %>% 
#   collect() 

#check that colnames match
cat(glue("Do all column headings match? {all(colnames(ind_swap) == colnames(ind_main))} \n"), " \n")

# Pre swap counts and checks
cat(glue("Number of individuals in db to be deleted: {sum(ind_main$study_id %in% swapout)} \n"), " \n")
cat(glue("Number of indiviuduals to be replaced: {sum(ind_swap$study_id %in% swapin)} \n"), " \n")
cat(glue("Did the swap change the count of individuals? {sum(ind_main$study_id %in% swapout) != sum(ind_swap$study_id %in% swapin)} \n"), " \n")
cat(glue("No replacement individuals currently in individual table? {(sum(ind_main$study_id %in% swapin))==0} \n"), " \n")

cat(glue("Number of new individuals to add to db: {sum(ind_swap$study_id %in% newstudies)} \n"), " \n")
cat(glue("No new individuals currently in db? {(sum(ind_main$study_id %in% newstudies))==0} \n"), " \n")

ind_final_pred <- nrow(ind_main) - sum(ind_main$study_id %in% swapout) + sum(ind_swap$study_id %in% swapin) + sum(ind_swap$study_id %in% newstudies)
cat(glue(" \n Final number of individuals should be: {cat(ind_final_pred)} \n"), " \n")

# sum(ind_loomis$study_id %in% swapout)
# sum(ind_loomis$study_id %in% swapin)

ind_out <- ind_main %>%  
  filter(study_id %notin% swapout) %>% # remove studies to be swapped
  bind_rows(ind_swap)   # bind new event table in

cat(glue("\n Final table length as expected? {nrow(ind_out)==ind_final_pred} \n"), " \n")

# write back to db
cat("Writing individual table back to db... \n")
dbWriteTable(dbmain, 'individual', ind_out, overwrite = T)


#---- Merge Event Tables ----#

cat("\n Merging event tables... \n")

# create lookup table to link 

evt_swap <- tbl(dbswap, 'event') %>% 
  left_join(ind_swap[,c("individual_id", "study_id")], copy = T) %>% # add study id to evt table
  mutate(study_id = as.character(study_id)) %>%
  collect()

# evt_loomis <- tbl(dbloomis, 'event_clean') %>% 
#   collect()

evt_main <- tbl(dbmain,'event') %>% # load main db event table
  left_join(ind_main[,c("individual_id", "study_id")], copy = T) %>% # add study id to evt table
  mutate(study_id = as.character(study_id)) %>% 
  collect()

#check that colnames match
cat(glue("Do all column headings match? {all(colnames(evt_swap) == colnames(evt_main))} \n"), " \n")

if(!all(colnames(evt_swap) == colnames(evt_main))){

  cat("Fixing columns... \n")
  #Fix columns...
colnames(evt_main)
colnames(evt_swap)

# Need to add 'population_d', 'stringency', 'county_fips', and 'census_tract' to swap table
evt_swap <- evt_swap %>% 
  mutate(population_d = NA,
         stringency= NA,
         county_fips = NA,
         census_tract = NA)
#check that colnames match
cat(glue("Do all column headings now match? {all(colnames(evt_swap) %in% colnames(evt_main))} \n"), " \n")
}

# Pre swap counts and checks
cat(glue("Number of events in db to be deleted: {sum(evt_main$study_id %in% swapout)} \n"), " \n")
cat(glue("Number of events to be replaced: {sum(evt_swap$study_id %in% swapin)} \n"), " \n")
cat(glue("Did the swap change the count of events? {sum(evt_main$study_id %in% swapout) != sum(evt_swap$study_id %in% swapin)} \n"), " \n")
cat(glue("No replacement individuals currently in event table? {(sum(evt_main$study_id %in% swapin))==0} \n"), " \n")
cat(glue("Number of new events to add to db: {sum(evt_swap$study_id %in% newstudies)} \n"), " \n")
cat(glue("No new events currently in db? {(sum(evt_main$study_id %in% newstudies))==0} \n"), " \n")

evt_final_pred <- nrow(evt_main) - sum(evt_main$study_id %in% swapout) + sum(evt_swap$study_id %in% swapin) + sum(evt_swap$study_id %in% newstudies)
cat(glue(" \nFinal number of events should be: {evt_final_pred} \n"), " \n")

evt_out <- evt_main %>%  
  filter(study_id %notin% swapout) %>% # remove studies to be swapped
  bind_rows(evt_swap)   # bind new event table in

cat(glue("\n Final table length as expected? {nrow(evt_out)==evt_final_pred} \n"), " \n")

# write back to db
cat("Writing event table back to db... \n")
dbWriteTable(dbmain, 'event', evt_out, overwrite = T)


cat("Disconnecting from db... \n")
dbDisconnect(dbmain)
dbDisconnect(dbswap)
# dbDisconnect(dbloomis)
# Export tables as csvs
# cat("Exporting data to csv...")
# write_csv(evt_main, file = "out/data_export/events.csv")
# write_csv(std_main, file = "out/data_export/study.csv")
# write_csv(ind_main, file = "out/data_export/individual.csv")


#---- Finalize script ----#

cat(glue('Script complete in {diffmin(t0)} minutes \n'), " \n")

sink()