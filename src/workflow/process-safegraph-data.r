#!/usr/bin/env Rscript --vanilla
# chmod 744 script_template.r #Use to make executable

# DESCRIPTION #
#
# This script reformats human mobility data for the COVID-19 Animal Movement Project
# See project documentation for details about anticipated directory structure.
#
# Major tasks fof this script:
#   * read in safegraph data
#   * reformat from wide to long format
#   * write out data by counties
#   * read in county files and combine into single file
#   * repeat for daily and hourly data

#
# This script implements the breezy philosophy: github.com/benscarlson/breezy

# ==== Breezy setup ====

#'
#Template
#Usage:
#script_template <taxa> <dat> <out> 
#script_template (-h | --help)
#Parameters:
#  dat: path to input csv file. 
#  out: path to output directory.
#Options:
#-h --help     Show this screen.
#-v --version     Show version.
#' -> doc

#---- Input Parameters ----#
if(interactive()) {
  rm(list=ls())
  library(here)
  
  .wd <- '~/Documents/Yale/projects/covid'
  .test <- TRUE
  rd <- here::here
  
  .datPF <- file.path(.wd,'data/safegraph/counties-dates-2-10-22/')
  .outPF <- file.path(.wd,'analysis/safegraph/counties-dates-2-10-22-reformatted/')
  
} else {
  library(docopt)
  library(rprojroot)
  
  #ag <- docopt(doc, version = '0.1\n')
  .wd <- '/gpfs/ysm/project/jetz/ryo3/projects/covid'
  .script <-  thisfile()
  #.test <- as.logical(ag$test)
  rd <- is_rstudio_project$make_fix_file(.script)
  
  #source(rd('src/funs/input_parse.r'))
  
  .datPF <- file.path(.wd,'data/safegraph/counties-dates-2-10-22/')
  .outPF <- file.path(.wd,'analysis/safegraph/counties-dates-2-10-22-reformatted/')
}

#---- Initialize Environment ----#

source(rd('src/startup.r'))

#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

# get list of files
files <- data.frame("name" = list.files(.datPF, pattern = "*.txt",recursive = TRUE)) %>%
  # separate file name into start date, county id, and resolution
  separate(name,into = c("start_date","county_id",NA,"resolution",NA), sep = "_", remove = FALSE) %>%
  separate(start_date, into = c("start_date", NA), sep = "/", remove = TRUE) %>%
  mutate("start_date" = lubridate::as_date(start_date))

# filter to daily data
files_daily <- files %>%
  filter(resolution == "day")

# filter to hourly data
files_hourly <- files %>%
  filter(resolution == "hour")

# get list of counties
counties <- unique(files$county_id)

# function for processing daily data
process_daily_data <- function(file_name){
  # read in file
  d <- fread(paste0(.datPF,file[,]$name), colClasses = "character") %>%
    # rename census block group column
    rename("cbg" = "V1") %>%
    # convert from "wide" to "long" format
    # columns are days of the week, values are device counts
    pivot_longer(!cbg, names_to = "day_of_week", values_to = "count", names_prefix = "V") %>%
    mutate("day_of_week" = as.numeric(day_of_week) - 2,
           "county_id" = rep(file[,]$county_id),
           "start_date" = rep(file[,]$start_date),
           "date" = start_date + day_of_week) %>% # find true date by adding from start date
    select(county_id,cbg,date,count)
  return(d)
}


for(i in 1:nrow(files_daily)){
  file <- files_daily[i,]
  d <- process_daily_data(paste0(.datPF,file[,]$name))
  fname <- str_sub(file$name, start = 12L, end = -5L)
  fwrite(d, paste0(.outPF,"daily-data/",fname,".csv"))
}






# DEPRECATED CODE

if (1 == 2){
  if(length(list.files(paste0(.outPF,"/daily-data/"))) == length(counties)){
    message("daily data already processed!")
  }else{
    for (j in 1:length(counties)){
      message("processing daily data...")
      
      # subset to files from single county
      county_files <- files_daily %>%
        filter(county_id == counties[j])
      
      d <- c()
      for (i in 1:nrow(county_files)){
        # start with first file
        file <- county_files[i,]
        # reformat data
        dd <- process_daily_data(file$name)
        d <- rbind(d,dd)
      }
      # write out file with all data from a single county
      fwrite(d, paste0(.outPF,"/daily-data/",counties[j],"_cbg_day_SUM.csv"))
    }
  } 
  
  if(file.exists(paste0(paste0(.outPF,"all_counties_cbg_day_SUM.csv")))){
    message("daily data already combined!")
  }else{
    message("combining daily data...")
    
    ## combine daily data into single file
    # combine county files into a single file
    reformatted_files_daily <- list.files(paste0(.outPF,"/daily-data/"), full.names = TRUE)
    
    # combine all data
    all_data_daily <- data.table::rbindlist(lapply(reformatted_files_daily, data.table::fread),use.names = TRUE)
    
    # write out as single file
    fwrite(all_data_daily, paste0(paste0(.outPF,"all_counties_cbg_day_SUM.csv")))
  }
  
  
  
  
  # create reformatted files for each county
  process_hourly_data <- function(file_name){
    # read in file
    d <- fread(paste0(.datPF,file[,]$name)) %>%
      # rename census block group column
      rename("cbg" = "V1") %>%
      # convert from "wide" to "long" format
      # columns are hours of the week, values are device counts
      pivot_longer(!cbg, names_to = "hour_of_week", values_to = "count", names_prefix = "V") %>%
      mutate("hour_of_week" = as.numeric(hour_of_week)-2,
             "day_of_week" = floor((as.numeric(hour_of_week))/24), #find day of week
             "hour_of_day" = hour_of_week - 24*day_of_week, #find hour of day
             "minutes_seconds" = rep(":00:00", nrow(.)), #create dummy minutes and seconds
             "county_id" = rep(file[,]$county_id),
             "start_date" = rep(file[,]$start_date),
             "date" = start_date + day_of_week) %>% # find true date by adding from start date
      unite("time", c("hour_of_day","minutes_seconds"),sep = "") %>% #concatenate time
      unite("date", c("date","time"), sep = " ") %>% #concatenate date and time
      select(county_id,cbg,date,count)
    return(d)
  }
  
  if(length(list.files(paste0(.outPF,"/hourly-data/"))) == length(counties)){
    message("hourly data already processed!")
  }else{
    message("processing hourl data...")
    
    # create reformatted files for each county
    for (j in 1:length(counties)){
      
      # subset to files from single county
      county_files <- files_hourly %>%
        filter(county_id == counties[j])
      
      d <- c()
      for (i in 1:nrow(county_files)){
        # start with first file
        file <- county_files[i,]
        # reformat data
        dd <- process_hourly_data(file$name)
        d <- rbind(d,dd)
      }
      # write out file with all data from a single county
      fwrite(d, paste0(.outPF,"/hourly-data/",counties[j],"_cbg_hour_SUM.csv"))
    }
  }
  
  if(file.exists(paste0(paste0(.outPF,"all_counties_cbg_hour_SUM.csv")))){
    message("hourly data already combined!")
  }else{
    message("combining hourly data...")
    
    ## combine hourly data into single file
    # combine county files into a single file
    reformatted_files_hourly <- list.files(paste0(.outPF,"/hourly-data/"), full.names = TRUE)
    
    # combine all data
    all_data_hourly <- data.table::rbindlist(lapply(reformatted_files_hourly, data.table::fread),use.names = TRUE)
    
    # write out as single file
    fwrite(all_data_hourly, paste0(paste0(.outPF,"all_counties_cbg_hour_SUM.csv")))
  }
  
}  