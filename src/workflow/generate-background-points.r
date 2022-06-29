if(interactive()) {
  # rm(list=ls())
  library(here)
  
  .wd <- '~/projects/covid-19_movement'
  .test <- TRUE
  # rd <- here::here
  .datPF <- file.path(.wd,'data/')
  .outPF <- file.path(.wd,'out/')
  
  .dbPF <- '/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db'
  
} else {
  library(docopt)
  library(rprojroot)
  
  .wd <- '/gpfs/ysm/project/jetz/ryo3/projects/covid'
  .script <-  thisfile()
  # rd <- is_rstudio_project$make_fix_file(.script)
  .datPF <- file.path(.wd,'data/')
  .outPF <- file.path(.wd,'out/')
  
  .dbPF <- '/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db'
  
}

message("start generating background points...")

source(file.path(.wd,'analysis/src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(data.table)
    library(lubridate)
    library(amt)
    library(DBI)
    library(RSQLite)
  }))



db <- dbConnect(RSQLite::SQLite(), .dbPF)



message("read in events and select random daily event...")
evt <- dbGetQuery(db,'SELECT event_id,individual_id,lon,lat,timestamp from event_clean') %>%
  collect() %>%
  mutate("date" = as_date(timestamp),
         "date_time" = as_datetime(timestamp),
         "year" = lubridate::year(date),
         "doy" = lubridate::yday(date)) %>%
  dplyr::filter(year >= 2019) %>% # lazy check on date filter
  dplyr::filter(year <= 2020) %>%
  dplyr::filter(doy < 170) %>%
  group_by(individual_id, date) %>%  
  sample_n(1) %>% # randomly select one event per day
  ungroup()



ids <- unique(evt$individual_id)
errors <- 0

for(i in 1:length(ids)){
  id <- ids[i]
  message(id)
  # filter to individual
  e <- evt %>%
    filter(individual_id == id) 
  
  # filter to individuals with greater than 10 records
  if (nrow(e) > 10){
    
    # generate track
    tr <- make_track(e, lon, lat, date_time, 
                     individual_id = individual_id)
    
    # CHECK TOLERANCE VALUE
    steps <- tr %>% 
      track_resample(rate = hours(24), tolerance = hours(22)) %>%
      steps_by_burst()
    
    # filter to individuals with steps (i.e. consistent data)
    if (nrow(steps) > 0){
      # check number of rows with NA turning angle
      test <- steps %>%
        filter(is.na(ta_))
      
      # filter to individuals with < 50% steps with NA turning angle
      if (nrow(test)/nrow(steps) < 0.5){
        
        # generate background points
        ssf <- steps %>%
          random_steps(n_control = 15) 
        
        fwrite(ssf, paste0(.outPF,"ssf-background-pts/individual-files/",id,".csv"))
        
        # write to log file
        write.table(c(ids[i], 1), glue("{.outPF}/sf-background-pts/bg-log.csv"), append = T, row.names = F, 
                    col.names = F, sep = ",")
        message(paste0(id,": worked"))
        
        test <- ssf %>%
          filter(abs(x2_) > 180) %>%
          filter(abs(y2_) > 90)
        
        if (nrow(test) > 0){
          message(paste0("errors (n): ", nrow(test)))
          errors <- errors + 1
        }
      } else{
        message(paste0(id,": filtered"))
      }
    } else{
      message(paste0(id,": filtered"))
    }
  } else{
    message(paste0(id,": filtered"))
  }
}

message(paste0("individuals with errors (n): ", errors))

dbDisconnect(db)

message("done!")