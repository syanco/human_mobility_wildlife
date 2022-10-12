################################################################################################################
#       COVID-19 Animal Movement                                                                               #
################################################################################################################

# This script generates individual n-dimensional hypervolumes (Lu et al. 2020).

# ==== Setup ====

'
Calculate multivariate niches across multiple individuals
Usage:

calc-niche-breadth.r <db> <out> <nc>
calc-niche-breadth.r (-h | --help)

Parameters:
  db: path to movement databse. 
  out: path to output directory.
  nc: number of cores for parallel processing
  
Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc


#---- Input Parameters ----#
if(interactive()) {
  # library(here)
  
  .wd <- getwd()
  
  .outPF <- file.path(.wd,'out/single_species_models/niche')
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
  
  .nc <- 2
  
} else {
  library(docopt)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  .nc <- ag$nc
  
}

#---- Initialize Environment ----#
# Start time
t0 <- Sys.time()

# Run startup
source(file.path(.wd,'analysis/src/startup.r'))

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
  }))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

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

# Look at climate data what is the percentage of NA. More than half of LST has NA
# prop.table(table(is.na(evt0$lst)))
# prop.table(table(is.na(evt0$tmax)))
# prop.table(table(is.na(evt0$tmin)))
# prop.table(table(is.na(evt0$ndvi)))
# prop.table(table(is.na(evt0$elev)))
# prop.table(table(is.na(evt0$dist2road)))


yearvec <- c("2019", "2020")

# Add empty columns study:

# log <- read_csv(glue("{.outPF}/niche_log.csv")) #####

registerDoMC(.nc)

# DEPRECATED - moved this process into the workflow.sh script
# # check if output csv exists and intit if not
# if(!file.exists(glue("{.outPF}/niche_determinant_anthropause_test.csv"))){
#   # Write an empty table:
#   determinant_template = data.frame()
#   columns= c("total","tmax","tmin","lst","ndvi", 'elev', 'dist2road', 'cor', 'week', 'individual',
#              'scientificname', 'studyid','year') 
#   determinant_template = data.frame(matrix(nrow = 0, ncol = length(columns))) 
#   # assign column names
#   colnames(determinant_template) = columns
#   
#   write.table(determinant_template, glue("{.outPF}/niche_determinant_anthropause_test.csv"), append = T, row.names = F, 
#               col.names = T, sep = ",")
# }

# DEPRECATED - moved this process into the workflow.sh script
# # check if log csv exists and intit if not
# if(!file.exists(glue("{.outPF}/niche_log_test.csv"))){
#   # Write the log file:
#   columns= c("studyid","individual","scientificname","year", 'status', 'week') 
#   logfile_template = data.frame(matrix(nrow = 0, ncol = length(columns))) 
#   # assign column names
#   colnames(logfile_template) = columns
#   
#   write.table(logfile_template, glue("{.outPF}/niche_log_test.csv"), append = T, row.names = F, 
#               col.names = T, sep = ",")
# }

foreach(j = 1:length(unique(ind)), .errorhandling = "pass", .inorder = F) %:%
  foreach(i = unique(yearvec), .errorhandling = "pass", .inorder = F) %dopar% {
# for(j in 1:length(unique(ind)) ){
#   # j <- 100
#   for(i in unique(yearvec)){
    # i <- 2019
    # print(paste0('Data for individual ', ind[j], ' year ', i))
    
    # Toggle `%do%` to `%dopar%` for HPC, %do% for local
    # foreach(j = 1:length(inds)) %do% {
    # foreach(j = 1:length(ind), .errorhandling = "pass", .inorder = F) %:%
    #   # j = 1
    #   foreach(i = 1:2, .errorhandling = "pass", .inorder = F) %dopar% { # do
    #     # i <- 2
    
    # check whether individual has been previously considered
    # if(log %>% 
    #    filter(ind_id == ind[j] & year == yearvec[i]) %>% 
    #    nrow() == 0){
    #   message( glue("Starting individual {ind[j]}, year {yearvec[i]}..."))  
    #   
    #---- Initialize database ----#
    
    #---- Perform analysis ----#
    message(glue("Gathering movement data for individual {ind[j]}, year {i}..."))
    
    scientificname <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(taxon_canonical_name)
    
    studyid <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(study_id)
    
    
    message("Filtering data and manipulating dates...")
    
    # print(i)
    
    # i = 2020
    evt_mod <- evt0 %>% 
      filter(individual_id == ind[j]) %>%
      dplyr::filter(yr == i) %>%
      mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %T"),
             week = week(timestamp),
             n_indiv_week_year = paste0(individual_id, '_' , week, '_' , yr),
             tmax_scale = scale(tmax),
             tmin_scale = scale(tmin),
             ndvi_scale = scale(ndvi),
             elev_scale = scale(elev),
             dist2road_scale = scale(dist2road)) %>% 
      arrange(timestamp)
    
    
    # tmp_sub = evt_mod %>% select("event_id","individual_id","timestamp","lon","lat", 'week', 'yr', 'n_indiv_week_year', 'tmax',
    #                              'tmin', 'ndvi', 'elev', 'dist2road' )
    
    
    #-- Fit Multivariate niches ####
    message("Estimating Multivariate niches")
    
    wk <- unique(evt_mod$week)
    # wk = wk[complete.cases(wk)]
    if(length(wk)==0){print(paste0('No data for year ', i, ' writing in logfile'))
      
      tmp_dummy_fail = data.frame(
        studyid = studyid, 
        individual  = ind[j],
        scientificname = scientificname,
        year = i,
        week = NA,
        status = 0)
      
      # logfile_template$studyid <- studyid
      write.table(tmp_dummy_fail, glue("{.outPF}/niche_log_test.csv"), append = T, row.names = F, 
                  col.names = F, sep = ",")
      
      
      } else{ # if no weeks in data
    
    
    # i <- 10
    for(w in min(wk):max(wk)){
      # w <- 26
      # Minimum of how many points per week? Maybe do biweekly? for sample size related? echoeing benas concerns? # Ask Scott ####
      # some weeks have no data NA values ind j == 4111
      # print(w)
      
      evt_tmp <- evt_mod %>% filter(week == w) %>%
        select(tmax_scale, tmin_scale, ndvi_scale, elev_scale, dist2road_scale, week, event_id, individual_id, 
               n_indiv_week_year) %>%
        drop_na(tmax, tmin, ndvi, elev, dist2road)
      
      # evt_tmp$n_indiv_week_year
      
      
      tryCatch({
        
        if(nrow(evt_tmp) > 0){      
          determinant <- MVNH_det(evt_tmp[,c('tmax_scale', 'tmin_scale', 'ndvi_scale', 'elev_scale', 'dist2road_scale')], log = F)

          determinant_df <- data.frame(as.list(determinant))
          determinant_df$week <- unique(w)
          determinant_df$individual <- unique(evt_tmp$individual_id)
          determinant_df$scientificname <- scientificname
          determinant_df$studyid <- studyid
          determinant_df$year <- i
          
          write.table(determinant_df, glue("{.outPF}/niche_determinant_anthropause_test_{Sys.Date}.csv"), append = T, row.names = F, 
                      col.names = F, sep = ",")
          
          
          
          tmp_dummy_success = data.frame(
            studyid = studyid, 
            individual  = ind[j],
            scientificname = scientificname,
            year = i,
            status = 1,
            week = w)
          
          # logfile_template$studyid <- studyid
          write.table(tmp_dummy_success, glue("{.outPF}/niche_log_test.csv"), append = T, row.names = F, 
                      col.names = F, sep = ",")
          } else {# if there's data
            
            tmp_dummy_success = data.frame(
              studyid = studyid, 
              individual  = ind[j],
              scientificname = scientificname,
              year = i,
              status = 0,
              week = w)
            
            # logfile_template$studyid <- studyid
            write.table(tmp_dummy_success, glue("{.outPF}/niche_log_test.csv"), append = T, row.names = F, 
                        col.names = F, sep = ",")
            
            print(paste0(unique(evt_tmp$n_indiv_week_year), ' has NA in niche determinant, moving to the next week'))
          } # else
            
      }, error = function(e){cat(glue("ERROR: unspecified error in fitting niche determinant for ind {ind[j]}, yr {i}", 
                                      "\n"))})
      
      
  
  } # for w in wks
      } # else (if wks > 0)
        #  } # fi end the check whether individual has been previously considered
} #i (end loop through years) : #j (end loop through individuals)


#---- Finalize script ----#

message(glue('Script complete in {diffmin(t0)} minutes'))
