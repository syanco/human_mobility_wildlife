library(here)
library(docopt)
library(rprojroot)

source(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/startup.r'))
source(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/funs/input_parse.r'))
#Source all files in the auto load funs directory
list.files(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/funs/auto'), full.names=T) %>% walk(source)
.dbPF <-'/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/Data/mosey_mod.db'
.wd <- '/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/'
.outPF <- paste0(.wd,'out/')

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
    require(MVNH)
    require(dplyr)
    require(doMPI)
  }))

#Source all files in the auto load funs directory
#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))
indtb <- tbl(db,'individual') %>%  # Load a tibble with all individual animals
  collect()

# Load the entire event table:
evt0 <- tbl(db, "event_clean")%>%
  collect()

evt0$timestamp_ymd = ymd_hms(evt0$timestamp)
evt0$jday = lubridate::yday(evt0$timestamp_ymd)
evt0$jday_year = paste0(evt0$jday, '-', evt0$yr)
evt0$individual_week_year = paste0(evt0$individual_id, '-',evt0$wk,'-', evt0$yr)

individual_week_year = unique(evt0$individual_week_year)
save(individual_week_year, file = '/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/indtb.Rdata')
# save(indiv_week_year, file = '/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/indtb.Rdata')
# save(data.frame(unique(evt0$individual_week_year)), file = '/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/indtb.Rdata')

