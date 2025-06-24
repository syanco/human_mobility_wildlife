#!/usr/bin/env Rscript --vanilla

# DESCRIPTION #
#
# This script reformats human mobility data for the COVID-19 Animal Movement Project
# See project documentation for details about anticipated directory structure.
#
# Major tasks fof this script:
#   * read in county files and at daily resolution
#   * read in safegraph data from the raw_data dir
#   * reformat from wide to long format
#   * write out data at daily resolution, by counties, to the processed_data dir

# This script implements the breezy philosophy: github.com/benscarlson/breezy

#---- Input Parameters ----#

library(data.table)

.wd <- getwd()
.datPF <- file.path(.wd,'raw_data/safegraph/counties-dates-2-10-22/')
.outPF <- file.path(.wd,'processed_data/safegraph/counties-dates-2-10-22-reformatted/')

#---- Initialize Environment ----#

source(file.path(.wd,'src/startup.r'))

# get list of county files
files <- data.frame("name" = list.files(.datPF, pattern = "*.txt",recursive = TRUE)) %>%
  # separate file name into start date, county id, and resolution
  separate(name,into = c("start_date","county_id",NA,"resolution",NA), sep = "_", remove = FALSE) %>%
  separate(start_date, into = c("start_date", NA), sep = "/", remove = TRUE) %>%
  mutate("start_date" = lubridate::as_date(start_date))

# filter to daily data
files_daily <- files %>%
  filter(resolution == "day")

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

# write out daily data with descriptive filenames
for(i in 1:nrow(files_daily)){
  file <- files_daily[i,]
  d <- process_daily_data(paste0(.datPF,file[,]$name))
  fname <- str_sub(file$name, start = 12L, end = -5L)
  fwrite(d, paste0(.outPF,"daily-data/",fname,".csv"))
}
