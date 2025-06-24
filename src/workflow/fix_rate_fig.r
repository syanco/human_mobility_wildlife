#!/usr/bin/env Rscript 
#
# This script generates fix rate for each ind/year combo, which can be used to
# determine the fix rater per study for supplementary figure 

# ==== Setup ====

'
Estimate fix rates for pre-segemented 
migratory periods across multiple individuals

Usage:
fix_rate_fig.r <db> <out> <nc>
fix_rate_fig.r (-h | --help)

Parameters:
  db: path to movement database. 
  out: path to output directory.
  nc: number of cores for parallel processing
  
Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  
  .wd <- "/home/julietcohen/repositories/human_mobility_wildlife"
  .dbPF <- file.path(.wd, 'processed_data/intermediate_db_copies/mosey_mod_clean-movement_complete.db')
  .outPF <- file.path(.wd,'out')
  
  .nc <- 1
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$db)
  .outPF <- makePath(ag$out)
  .nc <- ag$nc
  
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
    library(lubridate)
    library(move)
    library(doMC)
    library(foreach)
    library(glue)
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))
indtb <- tbl(db,'individual_clean') %>%
  mutate(taxon_canonical_name = case_when(
        study_id == 1442516400 ~ "Anser caerulescens",
        study_id == 1631574074 ~ "Ursus americanus",
        study_id == 1418296656 ~ "Numenius americanus",
        study_id == 2548691779 ~ "Odocoileus hemionus",
        study_id == 2575515057 ~ "Cervus canadensis",
        study_id == 1044238185 ~ "Alces alces",
        study_id == 474651680  ~ "Odocoileus virginianus",
        TRUE ~ taxon_canonical_name)) %>%
  mutate(taxon_canonical_name = case_when(
         taxon_canonical_name == "Chen rossii" ~ "Anser rossii",
         taxon_canonical_name == "Spilogale putorius" ~ "Spilogale interrupta",
         taxon_canonical_name == "Cervus elaphus" ~ "Cervus canadensis",
         taxon_canonical_name == "Accipiter gentilis" ~ "Astur atricapillus",
         taxon_canonical_name == "Chen caerulescens" ~ "Anser caerulescens",
         TRUE ~ taxon_canonical_name)) %>%
  # remove certain skunks per data owners request
  filter(individual_id != 3418130234, 
        individual_id != 3418130244, 
        individual_id != 3418130245, 
        individual_id != 3418130249) %>%
  collect() 

indtb <- indtb[!duplicated(indtb),]

message("Disconnecting from database...")
dbDisconnect(db)

ind <- indtb %>% 
  pull(individual_id)

#remove duplicated created due to merging of dbs
ind <- ind[!duplicated(ind)] 

yearvec <- c("2019", "2020")

registerDoMC(.nc)

message("Starting fix rate calcs")

foreach(j = 1:length(ind), .errorhandling = "pass", .inorder = F) %:%
  foreach(i = 1:2, .errorhandling = "pass", .inorder = F) %dopar% {

    message(glue("Starting individual {ind[j]}, year {yearvec[i]}..."))  
    
    invisible(assert_that(file.exists(.dbPF)))
    db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
    invisible(assert_that(length(dbListTables(db))>0))
    
    evt0 <- tbl(db, "event_clean")
    
    scientificname <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(taxon_canonical_name)
    
    studyid <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(study_id)
    
    evt_mod <- evt0 %>% 
      filter(individual_id == !!ind[j]) %>%
      collect() %>% 
      mutate(timestamp = base::as.POSIXct(timestamp, format = "%Y-%m-%d %T", tz = "UTC"),
             year = lubridate::year(timestamp)) %>%
      filter(year == !!yearvec[i]) %>% 
      distinct() %>% 
      # sort by timestamp
      arrange(timestamp)
    
    dbDisconnect(db)
    
    # check that the filtered dataframe has any data
    if (nrow(evt_mod) == 0) {

      message(glue("No paired data available for ind {ind[j]}, yr {yearvec[i]}; no fix rate available."))
      
      # move onto the next iteration of ind and year
      return(NULL)
      
    } else {

      tryCatch({
        # make minimal df for `move`    
        evt_tmp <- evt_mod %>% 
          select(lon, lat, timestamp, wk)
        
        evt_mv <- move(x=evt_tmp$lon, y=evt_tmp$lat, time=evt_tmp$timestamp,
                        proj=CRS("+proj=longlat"))

        # calculate time diff between each fix
        fixrates <- timeLag(x=evt_mv, units="mins")
        # calculate an excessive lag
        #fix_med_tripled <- (3*median(timeLag(x=evt_mv, units="mins")))
        # if max fix rate is greater than 3*median for this individual-year,
        # exclude it as we did when we calculated the dBBMM
        #fixrates_cleaned <- ifelse(fixrates > fix_med_tripled, NA, fixrates)
        fix_med <- median(fixrates, na.rm = T)
        fix_min <- min(fixrates, na.rm = T)
        fix_max <- max(fixrates, na.rm = T)

        fix_df <- data.frame("study_id" = studyid, 
                              "species" = scientificname,
                              "individual_id" = ind[j],
                              "year" = yearvec[i],
                              "fix_rate_med" = fix_med,
                              "fix_rate_min" = fix_min,
                              "fix_rate_max" = fix_max)

        write.table(fix_df, 
                glue("{.outPF}/fixrate_med_min_max.csv"), 
                append = T, 
                row.names = F, 
                col.names = F, 
                sep = ",")
      
      }, error = function(e){
        
        # catch the exact error reported from the move package
        error_msg <- conditionMessage(e)
        # log info for failure with error message
        message(glue("ERROR in estimating or saving fix rate for {scientificname} ind {ind[j]}, yr {yearvec[i]}:\n{error_msg}", 
                                      "\n"))
      })
    }
  }


library(tidyverse)
options(scipen = 999)

fixrate <- read_csv(glue("{.outPF}/fixrate_med_min_max.csv"))

fixrate_species <- fixrate %>%
  group_by(species) %>%
  # calc median fixrate of the medians across all studies 
  # within species grouping
  summarise(med_fixrate_minutes = median(fix_rate_med)) %>% 
  # convert minutes to hours
  mutate(med_fixrate_hours = med_fixrate_minutes/60)

write_csv(fixrate_species, glue("{.outPF}/fixrate_sp_median.csv"))

message(glue('Script complete in {diffmin(t0)} minutes'))
