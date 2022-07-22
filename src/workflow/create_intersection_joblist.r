if(interactive()) {
  rm(list=ls())
  library(here)
  
  .wd <- '/gpfs/ysm/project/jetz/ryo3/projects/covid'
  .test <- TRUE
  # rd <- here::here
  
  .dbPF <- '/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db'
  .datPF <- file.path(.wd,'src/workflow/')
  .outPF <- file.path(.wd,"src/workflow/")
  
} else {
  library(docopt)
  library(rprojroot)
  
  .wd <- getwd()
  .script <-  thisfile()
  # rd <- is_rstudio_project$make_fix_file(.script)
  
  .dbPF <- '/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db'
  .datPF <- file.path(.wd,'analysis/src/workflow/')
  .outPF <- file.path(.wd,"analysis/src/workflow/")
}

source(file.path(.wd,'analysis/src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(data.table)
  }))

invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

n_total <- 'select count(*) as num from event_clean' %>%
  dbGetQuery(db,.)


n <- 100
n_events <- floor(n_total$num/n)

start_ix <- seq(from = 1, to = n*n_events, by = n_events)
end_ix <- seq(from = n_events, to = n_total$num, by = n_events)


joblist <- data.frame("string" = rep(paste0(" module load R/4.1.0-foss-2020b; Rscript ",.datPF,"intersect-events-cbg.r "), times = n),
                      "arg1" = start_ix,
                      "arg2" = end_ix,
                      "arg3" = seq(from = 1, to = n, by = 1))

write.table(joblist,paste0(.outPF,"joblist.txt"),col.names = FALSE,row.names = FALSE,quote = FALSE)