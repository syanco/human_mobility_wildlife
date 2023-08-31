# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#
# Using detailed human activity and remote sensing data t assess wildlife responses to altered huan behavior during the COVID-19 pandemic
#
# diego.ellissoto@yale.edu
#
# The aim of this scirpt is to calculate niche breadth subsampling for 10, 20, 30 and 50 numbers of points and see how estimtes would dffer with changing sample sizes
#
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---



# NEED TO RENAME DETERMINANT INDIVIDUAL TO INTEGER ^$ NOT INTEGER SIMPLE -> EVT0 is integer 64 ! 


# Load Study ID into here: 
# require(bit64)
# load('/Users/diegoellis/projects/Anthropause/Supplementary_material/indtb.Rdata')
# # save(individual_week_year, file = '/Users/diegoellis/projects/Anthropause/Supplementary_material/indtb.Rdata')
# individual_week_year$individual = sub("_.*", "", individual_week_year$indiv_week_year)
# # head(sort(individual_week_year$individual))
# # head(sort(evt0$individual_id))
# 
# individual_week_year$individual = as.integer64(individual_week_year$individual)
#(individual_week_year[ unique(individual_week_year$individual) %in% evt0$individual_id,])

# class(evt0$individual_id)

# Remove weeks of the year that do not have cleaned data:
# indtb[!indtb$individual_id %in% evt0$individual_id]

# Keep only individuals with the data:
# (indtb[unique(indtb$individual_id) %in% unique(evt0$individual_id) ,])
# 
#CONDA: covid


#!/usr/bin/env Rscript 

# ==== Setup ====
# Calculate multivariate niches across multiple individuals Usage:
# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(here)
    library(docopt)
    library(rprojroot)
    library(move)
    library(lubridate)
    library(raster)
    library(doMC)
    library(foreach)
    require(dplyr)
    require(doMPI)
    require(stringr)
    require(tidyverse)
    require(MVNH)
  }))


# Run startup
if(interactive()){
  source(file.path('/Users/diegoellis/projects/Anthropause/src/startup.r'))  
  source(file.path('/Users/diegoellis/projects/Anthropause/src/funs/input_parse.r'))
  list.files(file.path('/Users/diegoellis/projects/Anthropause/src/funs/auto'), full.names=T) %>% walk(source)
  subsample_folders = list.files('/Users/diegoellis/projects/Anthropause/Supplementary_material/out/' , full.names=T, pattern = 'data_subset_anthropause')
  n_sample_folder = subsample_folders[4]
  # individual_week_year = read.csv('/Users/diegoellis/projects/Anthropause/Supplementary_material/individual_week_year.csv') %>% dplyr::select(-X)
  load('/Users/diegoellis/projects/Anthropause/Supplementary_material/indtb.Rdata')
  names(individual_week_year) <- 'indiv_week_year'
  .wd <- '//Users/diegoellis/projects/Anthropause/Supplementary_material/'
  .outPF <- paste0(.wd,'out/NVMH_subsampling/')
  .dbPF <- '/Users/diegoellis/Desktop/mosey_mod.db'
  if(!file.exists(.outPF)){dir.create(.outPF)}
  
}else{
  source(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/startup.r'))
  source(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/funs/input_parse.r'))
  #Source all files in the auto load funs directory
  list.files(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/funs/auto'), full.names=T) %>% walk(source)
  
  
  subsample_folders = list.files("/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/" , full.names=T, pattern = 'data_subset_anthropause')
  n_sample_folder = subsample_folders[4]
  # individual_week_year = read.csv('/gpfs/loomis/pi/jetz/de293/Anthropause/Code/Space_use_sensitivity/individual_week_year.csv') %>% select(-X)
  # individual_week_year = read.csv('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/individual_week_year.csv') %>% dplyr::select(-X)
  # names(individual_week_year) <- 'indiv_week_year'
  load('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/indtb.Rdata')
  .dbPF <-'/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/Data/mosey_mod.db'
  .wd <- '/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/'
  .outPF <- paste0(.wd,'out/NVMH_subsampling/')
  if(!file.exists(.outPF)){dir.create(.outPF)}
  }


#---- Initialize Environment ----#
# Start time
t0 <- Sys.time()
.nc <- 20
# .nc <- 1

if(!file.exists(n_sample_folder)){dir.create(n_sample_folder)}
# i <- 1

foreach(i = 1:nrow(individual_week_year), .errorhandling = 'pass', .inorder = F ) %dopar% {
  
  set.seed(1)
  
  message(paste0("Gathering movement data for individual, week, year ", individual_week_year[i,], ' with ', basename(n_sample_folder), ' number of points' ))
  
  tmp_file = individual_week_year$indiv_week_year[i]
  
  # Remove everything before first _
  individual_animal_id = gsub("\\_.*","",tmp_file) # escape everything after '_'
  year = gsub(".*_","",tmp_file)
  tmp_file = gsub(individual_animal_id, '', tmp_file) # remove individual animal id from the string
  tmp_file = substring(tmp_file, 2) # Remove _
  week = gsub("\\_.*","",tmp_file) # escape everything after '_'
  n_sample = gsub( 'data_subset_anthropause_', '', basename(n_sample_folder)) 
  
  # Differente el number del .csv al nombre del animal cuyos datos estan adentro !!!
  # Make a log table:
  # outlog <- data.frame(
  #   individual_animal_id = individual_animal_id,
  #   week = week,
  #   year = year,
  #   ind_animal_id_year_week_nsample = paste0(individual_animal_id, '-', year, '-', week, '-', n_sample),
  #   n_sample = n_sample,
  #   iteration = NA,
  #   total = NA,
  #   tmax = NA,
  #   tmin = NA,
  #   ndvi = NA,
  #   elev = NA,
  #   cor = NA
  # )


  evt_tmp = read.csv(
    paste0(n_sample_folder, '/data_subset_anthropause_', individual_week_year$indiv_week_year[i],'_', n_sample, '.csv')
  )%>%  # Scale  
    dplyr::mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %T"),
             week = week(timestamp),
             n_indiv_week_year = paste0(individual_id, '_' , week, '_' , yr),
             tmax_scale = scale(tmax),
             # tmin_scale = scale(tmin),
             ndvi_scale = scale(ndvi),
             elev_scale = scale(elev)) %>% 
    dplyr::arrange(timestamp)
    
     #-- Fit Multivariate niches ####
    message("Estimating 
            Multivariate niches")
 
  # Sample and calculate the niche breadth for these 100 iterations:
  # Sample and calculate the niche breadth for these 20 iterations:
  # for(j in 1:100){
    for(j in 1:20){
    # j<- 1
    # print(j)
    evt_tmp_sub = evt_tmp %>% dplyr::filter(iteration == j) %>% dplyr::arrange(timestamp)
    
    
    tryCatch({
          
          if(nrow(evt_tmp_sub) > 0){      
            # determinant <- MVNH_det(evt_tmp_sub[,c('tmax_scale', 'tmin_scale', 
            #                                    'ndvi_scale', 'elev_scale')], 
            #                         log = F)
            
            determinant <- MVNH_det(evt_tmp_sub[,c('tmax_scale', 
                                                   'ndvi_scale',
                                                   'elev_scale')], 
                                    log = F)
            
            determinant_df <- data.frame(as.list(determinant))
            determinant_df$week <- unique(evt_tmp_sub$week)
            determinant_df$individual <- unique(evt_tmp_sub$individual_id)
            # determinant_df$scientificname <- scientificname
            # determinant_df$studyid <- studyid
            determinant_df$year <- unique(evt_tmp_sub$yr)
            determinant_df$iteration = j
            determinant_df$n_sample = unique(evt_tmp_sub$n_samples)
            determinant_df$n_indiv_week_year  = unique(evt_tmp_sub$n_indiv_week_year)
            
            # Link to individual study ID 
            
            write.table(determinant_df, glue("{.outPF}niche_subsampling_50.csv"), append = T, 
                        row.names = F, col.names = F, sep = ",")
            
            # Write into log file
            tmp_dummy_success = data.frame(
              n_indiv_week_year  = unique(evt_tmp_sub$n_indiv_week_year),
              individual  = unique(evt_tmp_sub$individual_id),
              year = unique(evt_tmp_sub$yr),
              status = 1,
              week = unique(evt_tmp_sub$week)
              )

            # logfile_template$studyid <- studyid
            write.table(tmp_dummy_success, glue("{.outPF}niche_log_50.csv"), append = T, row.names = F, col.names = F, sep = ",")
          
          }else{
                       tmp_dummy_fail = data.frame(
              n_indiv_week_year  = unique(evt_tmp$n_indiv_week_year),
              individual  = unique(evt_tmp$individual_id),
              year = unique(evt_tmp$yr),
              status = 0,
              week = unique(evt_tmp$week)
              )
              
write.table(tmp_dummy_fail, glue("{.outPF}niche_log_50.csv"), append = T, row.names = F, col.names = F, sep = ",")
                        
              
          }
          
          
             }, error = function(e){cat(
          glue("ERROR: unspecified error in fitting niche determinant for ind {unique(evt_tmp_sub$individual_id)}")
          )
          })
          
          # If error add a fail into the log file

          # )})
            
            
          } # End o j
          } # End of foreach
    
      
  # write.table(determinant_template, glue("{.outPF}/niche_det_test_{n_sample}/niche_determinant_anthropause_{missing_indiv_week_year[i]}_{n_sample}.csv"), append = T, row.names = F, 
  #             col.names = T, sep = ",")
  
  
