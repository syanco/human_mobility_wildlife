# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#
# Subsampling of DBBMMs
#
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


#---- Input Parameters ----#z
suppressWarnings(
  suppressPackageStartupMessages({
    library(move)
    library(lubridate)
    library(raster)
    library(doMC)
    library(foreach)
    require(dplyr)
    require(doMPI)
    require(stringr)
    require(tidyverse)
  }))


if(interactive()){
  source(file.path('/Users/diegoellis/projects/Anthropause/src/startup.r'))  
  source(file.path('/Users/diegoellis/projects/Anthropause/src/funs/input_parse.r'))
  list.files(file.path('/Users/diegoellis/projects/Anthropause/src/funs/auto'), full.names=T) %>% walk(source)
  load('/Users/diegoellis/projects/Anthropause/Supplementary_material/indtb.Rdata')
  names(individual_week_year) <- 'indiv_week_year'
  .wd <- '/Users/diegoellis/projects/Anthropause/Supplementary_material/'
  .outPF <- paste0(.wd,'out/')
  .outPF <- '/Users/diegoellis/projects/Anthropause/Supplementary_material/out'
  .outPF <- '/Users/diegoellis/projects/Anthropause/out/Out_test_2/'
  # subsample_folders =  list.files('/Users/diegoellis/projects/Anthropause/out/Out_test_2/', full.names=T, pattern = 'data_subset_anthropause')
  subsample_folders =  list.files('/Users/diegoellis/projects/Anthropause/out/Out_test_2/data_subset_test_10', full.names=T)
  subsample_folders = '/Users/diegoellis/projects/Anthropause/out/Out_test_2/data_subset_test_10/'
  n_sample_folder = subsample_folders[2]
  # For local test
  if(Sys.getenv("USER")=='diegoellis'){
    n_sample = gsub( 'data_subset_test_', '', basename(n_sample_folder)) 
  }
  # n_sample_folder[1] = subsample_folders
  if(!file.exists(.outPF)){dir.create(.outPF)}
}else{
  source(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/startup.r'))
  source(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/funs/input_parse.r'))
  list.files(file.path('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/analysis/src/funs/auto'), full.names=T) %>% walk(source)
  subsample_folders = list.files("/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/" , full.names=T, pattern = 'data_subset_anthropause')
  n_sample_folder = subsample_folders[2]
  load('/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/out/indtb.Rdata')
  names(individual_week_year) <- 'indiv_week_year'
  .dbPF <-'/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/Data/mosey_mod.db'
  .wd <- '/gpfs/gibbs/pi/jetz/from_loomis/de293/Anthropause/Anthropause_2023/'
  .outPF <- paste0(.wd,'out/')
  if(!file.exists(.outPF)){dir.create(.outPF)}
  n_sample = gsub( 'data_subset_anthropause_', '', basename(n_sample_folder)) 
}


#---- Initialize Environment ----#
# Start time
t0 <- Sys.time()
.nc <- 15
# .nc <- 1
if(!file.exists(n_sample_folder)){dir.create(n_sample_folder)}


foreach(i = 1:nrow(individual_week_year), .errorhandling = 'pass', .inorder = F ) %dopar% {
  # i=1
  set.seed(1)
  message(paste0("Gathering movement data for individual, week, year ", individual_week_year[1, i], ' with ', basename(n_sample_folder), ' number of points' ))
  
  tmp_file = individual_week_year$indiv_week_year[i]
  
  # Remove everything before first _
  individual_animal_id = gsub("\\_.*","",tmp_file) # escape everything after '_'
  year = gsub(".*_","",tmp_file)
  tmp_file = gsub(individual_animal_id, '', tmp_file) # remove individual animal id from the string
  tmp_file = substring(tmp_file, 2) # Remove _
  week = gsub("\\_.*","",tmp_file) # escape everything after '_'
  
  
  
  # Differente el number del .csv al nombre del animal cuyos datos estan adentro !!!
  # Make a log table:
  outlog <- data.frame(
    individual_animal_id = individual_animal_id,
    n_sample = n_sample,
    week = week,
    year = year,
    ind_animal_id_year_week_nsample = paste0(individual_animal_id, '-', year, '-', week, '-', n_sample),
    iteration = NA,
    dbmm_size = NA
  )
  
  if( ! file.exists(glue("{.outPF}/DBBMM_subsampling/dbbmm_20_log.csv"))  ){
  write.table(
    data.frame(
      individual_animal_id = NA,
      n_sample = NA,
      week = NA,
      year = NA,
      ind_animal_id_year_week_nsample = NA,
      iteration = NA,
      dbmm_size = NA
    )
    
    , glue("{.outPF}/DBBMM_subsampling/dbbmm_20_log.csv"), append = T, row.names = F,col.names = T, sep = ",") # Log file
  }
  
  if(Sys.getenv("USER")=='diegoellis'){
    evt_tmp = read.csv(
      paste0(n_sample_folder, 'data_subset_anthropause_',individual_week_year$indiv_week_year[i], '_',n_sample,'.csv') )%>%
      dplyr::select(event_id, individual_id,
                    timestamp, lon, lat, week, yr,
                    n_indiv_week_year, iteration,
                    n_samples,n_indiv_week_year_iteration)
    
  }else{
    # For HPC: 
    evt_tmp = read.csv(
      paste0(n_sample_folder,'/data_subset_anthropause_' ,individual_week_year$indiv_week_year[i],'_', n_sample, '.csv') )%>%
      dplyr::select(event_id, individual_id,
                    timestamp, lon, lat, week, yr,
                    n_indiv_week_year, iteration,
                    n_samples,n_indiv_week_year_iteration)
  }
  
  # Loop one hundred times to make dbmm n_sample
  # Sample and calculate the space use for these 100 iterations:
  for(j in 1:20){
    # j<- 1
    # print(j)
    evt_tmp_sub = evt_tmp %>% dplyr::filter(iteration == j) %>% dplyr::arrange(timestamp)
    
    evt_mv <- move(x=evt_tmp_sub$lon,
                   y=evt_tmp_sub$lat,
                   time=as.POSIXct(evt_tmp_sub$timestamp, format="%Y-%m-%d %H:%M:%OS", tz="UTC"), # trt = evt_tmp_sub$trt,
                   proj=CRS("+proj=longlat"),
                   data = evt_tmp_sub)
    
    evt_mv$week <- week(evt_mv$timestamp)
    burstid <- factor(evt_mv$week[1:(n.locs(evt_mv)-1)])
    #id "intended fix rate"
    fixmed <- median(timeLag(x=evt_mv, units="mins"))
    evt_burst <- burst(evt_mv, burstid)
    evt_mv_t <- spTransform(evt_burst, center = T)
    
    # Margin and window size need to be odd:
    margin_size = (as.numeric(n_sample) / 2)  - 1 # Cant be just half the sample size because margin to window ratio will not be appropiate
    
    if((margin_size %% 2) == 0){
      margin_size_tmp <- margin_size - 1
    }
    
    tryCatch({
      dbb_var <- brownian.motion.variance.dyn(
        object = evt_mv_t,
        location.error=5,
        margin = margin_size_tmp, # Margin divided -1 and round to closest smaller uneven integer
        window.size = (n.locs(evt_mv_t) - 1 )
      )
      
      # remove any segments with gaps > 3x the intended fix rate
      dbb_var@interest[timeLag(evt_mv_t,"mins")>(fixmed*3)] <- FALSE
      
      # 10x the raster size of the underlying data 
      dbbm <- brownian.bridge.dyn(
        dbb_var, 
        time.step = (fixmed/15),
        location.error = rep(5, n.locs(evt_mv_t)),
        ext = 10,
        margin = margin_size_tmp,
        window.size = (n.locs(evt_mv_t) - 1 )
      )
      
      tmp_out <- list("dBBMM Variance" = dbb_var,
                      "dBBMM Object" = dbbm,
                      "events" = evt_tmp_sub
      )

      
      r <- tmp_out$`dBBMM Object`
      rb <- UDStack(r)
      UDr <- getVolumeUD(rb)
      ud95 <- UDr<=.95
      a <- sum(values(ud95))*res(ud95)[1]*res(ud95)[1]
      
      
      outlog$dbmm_size <- a
      outlog$iteration <- j
      # Make entry in log file
      message(glue("Writing output for individual {individual_animal_id} to file..."))
        write.table(outlog, glue("{.outPF}/DBBMM_subsampling/dbbmm_20_log.csv"), append = T, row.names = F,col.names = F, sep = ",") # Log file
      
      
    }, error = function(e){cat(glue("ERROR: unspecified error in fitting dBBMM for ind {individual_animal_id}", 
                                    "\n"))})
    
    
    # message(glue("Writing output for individual {gsub('.csv','',gsub('data_subset_anthropause_','',basename(n_sample_folder[i])))} to file..."))

    # if(exists("outlog")){rm(outlog)}
    if(exists("tmp_out")){rm(tmp_out)}
    if(exists("dbbm")){rm(dbbm)}
    
  } # End of j iterations subsampling
  
} # End of loop individuals
