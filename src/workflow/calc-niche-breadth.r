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
  
  .outPF <- file.path(.wd,'out/niche_determinant_anthropause.csv')
  .dbPF <- file.path(.wd,'processed_data/mosey_mod_2023.db')
  
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
indtb <- tbl(db,'individual_final') %>%  # Load a tibble with all individual animals
  collect()

# Load the entire event table:
evt0 <- tbl(db, "event_final")%>%
  collect()

message("Disconnecting from databse...")
dbDisconnect(db)

ind <- indtb %>% 
  pull(individual_id)

ind <- ind[ind  %in% unique(evt0$individual_id)]

yearvec <- c("2019", "2020")

registerDoMC(.nc)

foreach(j = 1:length(unique(ind)), .errorhandling = "pass", .inorder = F) %:%
  foreach(i = unique(yearvec), .errorhandling = "pass", .inorder = F) %dopar% {
    
    
    #---- Perform analysis ----#
    message(glue("Gathering movement data for individual {ind[j]}, year {i}..."))
    
    scientificname <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(taxon_canonical_name)
    
    studyid <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(study_id)
    
    
    message("Filtering and scaling data...")
    
    # print(i)
    
    # i = 2020
    evt_mod <- evt0 %>% 
      filter(individual_id == ind[j]) %>%
      dplyr::filter(yr == i) %>%
      mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %T"),
             week = week(timestamp),
             n_indiv_week_year = paste0(individual_id, '_' , week, '_' , yr),
             tmax_scale = scale(tmax),
             ndvi_scale = scale(ndvi),
             elev_scale = scale(elev)
             ) %>% 
      arrange(timestamp)
    
    
    #-- Fit Multivariate niches ####
    message("Estimating Multivariate niches")
    
    wk <- unique(evt_mod$week)
    if(length(wk)==0){print(paste0('No data for year ', i, 
                                   ' writing in logfile'))
      
      tmp_dummy_fail = data.frame(
        studyid = studyid, 
        individual  = ind[j],
        scientificname = scientificname,
        year = i,
        week = NA,
        status = 0)
      
      write.table(tmp_dummy_fail, glue("./out/niche_log.csv"), append = T, 
                  row.names = F, col.names = F, sep = ",")
      
    } else { # if no weeks in data
      
      # i <- 10
      for(w in min(wk):max(wk)){
        # w <- 26
        
        evt_tmp <- evt_mod %>% 
          filter(week == w) %>%
          select(tmax_scale, ndvi_scale, elev_scale, week, event_id,
                 individual_id, n_indiv_week_year) %>%
          drop_na(tmax_scale, ndvi_scale, elev_scale)
        
        
        tryCatch({
          
          if(nrow(evt_tmp) > 0){      
            determinant <- MVNH_det(evt_tmp[,c('tmax_scale', 'ndvi_scale', 'elev_scale')], 
                                    log = F)
            
            determinant_df <- data.frame(as.list(determinant))
            determinant_df$week <- unique(w)
            determinant_df$individual <- unique(evt_tmp$individual_id)
            determinant_df$scientificname <- scientificname
            determinant_df$studyid <- studyid
            determinant_df$year <- i
            
            write.table(determinant_df, glue("{.outPF}"), append = T, 
                        row.names = F, col.names = F, sep = ",")
            
            tmp_dummy_success = data.frame(
              studyid = studyid, 
              individual  = ind[j],
              scientificname = scientificname,
              year = i,
              status = 1,
              week = w)
            
            # logfile_template$studyid <- studyid
            write.table(tmp_dummy_success, glue("./out/niche_log.csv"), 
                        append = T, row.names = F, 
                        col.names = F, sep = ",")
          } else {# if there's no data
            
            tmp_dummy_fail = data.frame(
              studyid = studyid, 
              individual  = ind[j],
              scientificname = scientificname,
              year = i,
              status = 0,
              week = w)
            
            # logfile_template$studyid <- studyid
            write.table(tmp_dummy_fail, glue("./out/niche_log.csv"), 
                        append = T, row.names = F, 
                        col.names = F, sep = ",")
            
            print(paste0(unique(evt_tmp$n_indiv_week_year), 
                         ' has NA in niche determinant, moving to the next week'))
          } # else
          
        }, error = function(e){cat(
          glue("ERROR: unspecified error in fitting niche determinant for ind {ind[j]}, yr {i}", 
                                        "\n"))})
        
      } # for w in wks
    } # else (if wks > 0)
    #  } # fi end the check whether individual has been previously considered
  } #i (end loop through years) : #j (end loop through individuals)


#---- Finalize script ----#

message(glue('Script complete in {diffmin(t0)} minutes'))
