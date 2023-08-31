# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#
# Subsampling of data: 
#
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#!/usr/bin/env Rscript 
#
#CONDA: covid

#---- Input Parameters ----#

library(here)
library(docopt)
library(rprojroot)


rd <- here::here
# .wd <- getwd()

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

#---- Initialize Environment ----#
t0 <- Sys.time()
# Run startup
if(interactive()){
  source(file.path('/Users/diegoellis/projects/Anthropause/src/startup.r'))  
  source(file.path('/Users/diegoellis/projects/Anthropause/src/funs/input_parse.r'))
  list.files(file.path('/Users/diegoellis/projects/Anthropause/src/funs/auto'), full.names=T) %>% walk(source)
  .dbPF <- '/Users/diegoellis/Desktop/mosey_mod.db'
  .wd <- '/Users/diegoellis/projects/Anthropause/Supplementary_material/'
  .outPF <- paste0(.wd,'out/')
}else{
  source(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/startup.r'))
  source(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/funs/input_parse.r'))
  #Source all files in the auto load funs directory
  list.files(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/funs/auto'), full.names=T) %>% walk(source)
  .dbPF <-'/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/Data/mosey_mod.db'
  .wd <- '/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/'
  .outPF <- paste0(.wd,'out/')
}


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

message("Disconnecting from databse...")
dbDisconnect(db)

ind <- indtb %>% 
  pull(individual_id)

ind = ind[ind  %in% unique(evt0$individual_id)]

# Create index
evt0$week <- week(evt0$timestamp)
evt0$n_indiv_week_year <- paste0(evt0$individual_id,'_' ,evt0$week, '_', evt0$yr)

# evt0 = read.csv('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/subset_evt0.csv')
# keep only the first 100 unique individual_id values
# example_tibble <- evt0 %>% 
#   group_by(individual_id) %>% 
#   filter(row_number() <= 20) %>% 
#   ungroup()
# load('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/indtb.Rdata')

# Look at the missing weeks:
unique_indiv_week_year <- unique(evt0$n_indiv_week_year)
indiv_week_year  = unique(evt0$n_indiv_week_year)
#---- Perform analysis ----#

# n_sample = 50
# n_sample = 30
# n_sample = 20
 n_sample = 10


# .nc <- 3
.nc <- 20
message(glue("Starting {.nc} cores"))
registerDoMC(.nc)

# i <- 2

# Check indiv_week_year 
foreach(i = 1:length(indiv_week_year), .errorhandling = 'pass',
        .inorder = F ) %dopar% {
          
          # If .csv file of this specific individual animal week year exists, next:
          if(
            file.exists(
              glue("{.outPF}data_subset_anthropause_{n_sample}/data_subset_anthropause_{indiv_week_year[i]}_{n_sample}.csv")
            )
          ){print('File exists going to next');next}
          
          
          
          set.seed(1)
          message(paste0("Gathering movement data for individual, week, year ", indiv_week_year[i]))
          
          scientificname <- indtb %>% 
            filter(individual_id == sub("_.*", "", indiv_week_year[i]) ) %>%  # individual_week_year
            pull(taxon_canonical_name)
          
          studyid <- indtb %>% 
            filter(individual_id == sub("_.*", "", indiv_week_year[i]) ) %>% 
            pull(study_id)
          
          if(!file.exists(glue("{.outPF}data_subset_anthropause_{n_sample}/"))){
            dir.create(glue("{.outPF}data_subset_anthropause_{n_sample}/"))
          }

          # No need to print this?          
          if(
            !file.exists(
             glue("{.outPF}data_subset_anthropause_{n_sample}/data_subset_anthropause_{indiv_week_year[i]}_{n_sample}.csv")
            )
          ){ # for each unique week of the year of each individual
          #   
            # Write an empty table:
            data_subset_template = data.frame()
            columns= c("event_id","individual_id","timestamp",
                       "lon","lat", 'week', 'yr',
                       'n_indiv_week_year', 'tmax','ndvi', 'elev',
                       'iteration', 'n_indiv_week_year_iteration', 'n_samples')
            
            data_subset_template = data.frame(matrix(nrow = 0, ncol = length(columns)))
            # assign column names
            colnames(data_subset_template) = columns
            
            # 
            
            
            # 
            write.table(data_subset_template,
                        file = glue("{.outPF}data_subset_anthropause_{n_sample}/data_subset_anthropause_{indiv_week_year[i]}_{n_sample}.csv"),
                        append = T, row.names = F,
                        col.names = T, sep = ",")
          
            # glue("{.outPF}data_subset_anthropause_{n_sample}/data_subset_anthropause_{indiv_week_year[i]}_{n_sample}.csv")
            }
          
          
          # Iterate 100 times:
          for(j in 1:20){
            print(j)
            #j = 1
            
            tryCatch({
              
              tmp_sub = evt0 %>% 
                dplyr::filter(n_indiv_week_year == indiv_week_year[i]) %>% 
                drop_na(tmax, ndvi, elev) %>%
                dplyr::select(event_id, individual_id, timestamp, lon, lat,
                              week, yr, n_indiv_week_year, tmax,
                              ndvi, elev) %>% 
                sample_n(size = n_sample,replace = FALSE)
              
              tmp_sub$iteration = j
              tmp_sub$n_indiv_week_year_iteration <- paste0(tmp_sub$n_indiv_week_year, '_n_', j)
              tmp_sub$n_samples = n_sample
              
              tmp_sub = tmp_sub %>% dplyr::select("event_id","individual_id","timestamp",
                                                  "lon","lat", 'week', 'yr',
                                                  'n_indiv_week_year', 'tmax',
                                                  'ndvi', 'elev', 'iteration',
                                                  'n_indiv_week_year_iteration', 'n_samples' )
              
              
              # Save subset of data
              write.table(tmp_sub, glue("{.outPF}data_subset_anthropause_{n_sample}/data_subset_anthropause_{indiv_week_year[i]}_{n_sample}.csv"),
                          append = T, row.names = F,
                          col.names = F, sep = ",")
              
            }, error = function(e){paste0(glue("ERROR: data_subset_anthropause_{indiv_week_year[i]}_{n_sample}"), ' has less than ', n_sample, ' points. Moving to next week'
            )})
            
            
          } # End j
          
        } # End i

