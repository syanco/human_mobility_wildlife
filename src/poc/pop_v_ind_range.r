###########################
#                         #
# Pop vs Ind Niche Dissim #
#                         #
###########################

# Scratch script to test out ideas of comapring population niche dissim with individual.

# LIBRARIES

library(MVNH)
library(jsonlite)
library(httr)
library(RSQLite)
library(DBI)
library(tidyverse)
library(lubridate)
library(glue)

# FUNCTIONS

get_scenarios <- function () {
  resp <- httr::GET(paste0(BASE_URL, '/list/scenarios'))
  jsonlite::flatten(jsonlite::fromJSON(httr::content(resp, "text", encoding = enc), simplifyVector = T))
}

get_species_scenario <- function (sciname, product, variable, s_buff, t_buff, limit=50000) {
  resp <- httr::POST(paste0(BASE_URL, '/species/metrics/scatter'), body = list(
    mode = 'temporal',
    scientificname = sciname,
    limit= limit,
    yaxis = list(
      product = product,
      variable = variable,
      temporal = t_buff, 
      spatial = s_buff
    )
  ), encode = "json",)
  jsonlite::flatten(jsonlite::fromJSON(httr::content(resp, "text", encoding = enc))$rows)
}

`%notin%` <- Negate(`%in%`)

# INITS

BASE_URL <- 'https://stoat-otf.azurewebsites.net'
enc = "UTF-8"
.dbPF <- file.path('data/anno_move.db')
.anno <- "anno_join_2021-09-08"
.ctfs <- file.path("ctfs/crane_segmentation")


# Pull STOAT Pre-Annotations
#
# View possible annotations
get_scenarios()

get_species_scenario('Grus grus', 'landsat8', 'evi', 250, 16)

#(while STOAT is doen...)
# Pull in Pop scale niche data
pop <- read.csv("analysis/cranes/Anthropoides virgo_landsat8_evi_500_16.csv")


# Read in cranes

#get segmentation filenames
segs <- list.files(.ctfs, pattern = "*.csv")

db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)

# Get target inds from one sp
inds <- tbl(db, "individual") %>%
  filter(taxon_canonical_name == "Anthropoides virgo") %>% 
  pull(individual_id)

# create df of annotations
anno0 <- tbl(db, .anno) %>% 
  filter(individual_id %in% inds) %>%     
  filter(is.na(`ground_speed`) | `ground_speed`<10) %>% 
  filter(is.na(gps_hdop) | gps_hdop < 5) %>% 
  collect()


# Annotate seasons from the anno files
# j <- 15 for tests
for(j in 1:length(ind)){
  message(glue("Starting individual {inds[j]}..."))
  
  message("Filtering data and manipulating dates...")
  
  # extract ind
  evt_mod <- anno0 %>% 
    filter(individual_id == inds[j]) %>% 
    mutate(date = date(timestamp),
           yr = year(date))
  
  
  #-- Load Segmentation Info
  message("Getting seasonal segmentations...")
  
  # grab the segmentation file by matching filename to ind
  if(inds[j] %in% str_extract(segs, "[0-9]+")){ #if the file exists...
    seg_temp <- read_csv(file.path( # ...read it in...
      .ctfs,segs[which(str_extract(segs, "[0-9]+") == inds[j])[1]] # ...mathcin individual_id
    )) %>%
      arrange(Date) %>% # sort by date
      #remove stopovers
      filter(Status %in% c("Start Fall", "End Fall", "Start Spring", "End Spring")) %>% 
      mutate(yr = year(Date)) # add year col
    
  } else {
    next
  }
  
  
  # add a variable for the stationary code
  seg_full <- seg_temp %>% 
    mutate(res_period = case_when(Status == "End Fall" ~ "Winter",
                                  Status == "Start Spring" ~ "Winter",
                                  Status == "End Spring" ~ "Summer",
                                  Status == "Start Fall" ~ "Summer"),
           mig_period = case_when(Status == "Start Fall" ~ "Fall",
                                  Status == "End Fall" ~ "Fall",
                                  Status == "Start Spring" ~ "Spring",
                                  Status == "End Spring" ~ "Spring")) %>% 
    arrange(Date)
  
  # TODO:  This process will still allow segmentation files with the correct
  #   sequence of stop-start, but across years, so encompassing 5 seasons
  #   need to build catch for that Grus grus (sp=3) ind=7 has this situation
  
  # Winter
  seg_res <- seg_full %>%
    # {if(nrow(.) < 2){.}else{
    #   # must start with season "openning" transition
    #   {if(.$Status[1] != "End Fall"){slice(., -1)}else{.}} %>% 
    #     # and end with season closing transition
    #     {if(.$Status[nrow(.)] != "Start Spring"){slice(., -nrow(.))}else{.}} %>% 
        # create repeating id columnto assist the pivot
        {if(nrow(.) > 1){
          mutate(., id = rep(seq(1, nrow(.)/2, by = 1), each = 2))}else{
            .
          }} %>%
        # remove unneeded columns
        select(-c(X1)) %>% 
        # pivot to wide format (1 row per season)
        pivot_wider(values_from = Date, names_from = Status) 
seg_res
  
  
  
  # TODO: pickup here
    
    
    #check for segmentation file data - set season to NA if segmentation is empty
    if(nrow(seg_temp) < 1){
      message("Insufficient segmentation data, moving to next individual!")
      evt_out <- evt_yr %>% 
        mutate(season = NA)
    }else{
      
      # Annotate the migratory periods first
      
      # Fall
      for(y in unique(seg_res$yr)){
        seg_res %>% 
          filter(yr == y)
      }
      
    }
    
    
  }
  
  
  
  evt_mod_season <- evt_mod %>% 
    ,
  season = case_when(date >= seg_temp$Date[seg_temp$Status == "Start Fall"] & date <= seg_temp$Date[seg_temp$Status == "End Fall"] ~ "Fall")),
date > seg_temp$Date[seg_temp$Status == "End Fall"] & date < seg_temp$Date[seg_temp$Status == "Start Spring"] ~ "Summer",
date >= seg_temp$Date[seg_temp$Status == "Start Spring"] & date <= seg_temp$Date[seg_temp$Status == "End Spring"] ~ "Spring"))


# add a variable for the stationary code
seg_full <- seg_temp %>% 
  mutate(res_period = case_when(Status == "End Fall" ~ "Winter",
                                Status == "Start Spring" ~ "Winter",
                                Status == "End Spring" ~ "Summer",
                                Status == "Start Fall" ~ "Summer"),
         mig_period = case_when(Status == "Start Fall" ~ "Fall",
                                Status == "End Fall" ~ "Fall",
                                Status == "Start Spring" ~ "Spring",
                                Status == "End Spring" ~ "Spring")) %>% 
  arrange(Date)



#-- Winter Niches
message("Gathering individual stats for winter period(s)...")

# init empty list for results
wint_out <- list()

if(nrow(seg_wint) > 0 & "run" %in% colnames(seg_wint)){
  for(i in 1:nrow(seg_wint)){
    
    # extract season of interest
    evt_tmp <- evt_mod %>% 
      filter(timestamp > seg_wint$`End Fall`[i] & timestamp < seg_wint$`Start Spring`[i])
    
    # Move to next season in the loop
    if(nrow(evt_tmp) == 0){
      message("No suitable records during season, moving to next...")
      next 
    }  
    
    # write mean and var to tmp
    tmp_out <- data.frame(
      season = "Winter",
      ind = as.character(ind[j]),
      mean = mean(na.omit(evt_tmp$`value_derived:evi`)),
      var = var(na.omit(evt_tmp$`value_derived:evi`)),
      n = nrow(evt_tmp)
    )
    wint_out[[i]] <- tmp_out
  } #i
} # fi      

#-- Summer Niches
message("Gathering individual stats for summer period(s)...")

# init empty list for results
summ_out <- list()

if(nrow(seg_summ) > 0 & "run" %in% colnames(seg_summ)){
  for(i in 1:nrow(seg_summ)){
    
    #format data as `telemetry` object for ctmm
    evt_tmp <- evt_mod %>% 
      filter(timestamp > seg_summ$`End Spring`[i] & timestamp < seg_summ$`Start Fall`[i])
    
    # Move to next season in the loop
    if(nrow(evt_tmp) == 0){
      message("No suitable records during season, moving to next...")
      next 
    }  
    
    # write mean and var to tmp
    tmp_out <- data.frame(
      season = "Summer",
      ind = as.character(ind[j]),
      mean = mean(na.omit(evt_tmp$`value_derived:evi`)),
      var =var(na.omit(evt_tmp$`value_derived:evi`)),
      n = nrow(evt_tmp)
    )
    summ_out[[i]] <- tmp_out
  } #i
} # fi   

#-- Spring Niches
message("Gathering individual stats for spring period(s)...")

# init empty list for results
spring_out <- list()

if(nrow(seg_spring) > 0 & "run" %in% colnames(seg_spring)){
  for(i in 1:nrow(seg_spring)){
    
    #format data as `telemetry` object for ctmm
    evt_tmp <- evt_mod %>% 
      filter(timestamp > seg_spring$`Start Spring`[i] & timestamp < seg_spring$`End Spring`[i])
    
    # Move to next season in the loop
    if(nrow(evt_tmp) == 0){
      message("No suitable records during season, moving to next...")
      next 
    }  
    
    # write mean and var to tmp
    tmp_out <- data.frame(
      season = "Spring",
      ind = as.character(ind[j]),
      mean = mean(na.omit(evt_tmp$`value_derived:evi`)),
      var =var(na.omit(evt_tmp$`value_derived:evi`)),
      n = nrow(evt_tmp)
    )
    spring_out[[i]] <- tmp_out
  } #i
} # fi   

#-- Fall Niches
message("Gathering individual stats for fall period(s)...")

# init empty list for results
fall_out <- list()

if(nrow(seg_fall) > 0 & "run" %in% colnames(seg_fall)){
  for(i in 1:nrow(seg_fall)){
    
    #format data as `telemetry` object for ctmm
    evt_tmp <- evt_mod %>% 
      filter(timestamp > seg_fall$`Start Fall`[i] & timestamp < seg_fall$`End Fall`[i])
    
    # Move to next season in the loop
    if(nrow(evt_tmp) == 0){
      message("No suitable records during season, moving to next...")
      next 
    }  
    
    # write mean and var to tmp
    tmp_out <- data.frame(
      season = "Fall",
      ind = as.character(ind[j]),
      mean = mean(na.omit(evt_tmp$`value_derived:evi`)),
      var =var(na.omit(evt_tmp$`value_derived:evi`)),
      n = nrow(evt_tmp)
    )
    fall_out[[i]] <- tmp_out
  } #i
} # fi   

#-- Collate outputs
message("Collating output...")

ind_ls <- list(
  wint = if(length(wint_out) > 0){do.call("rbind", wint_out)}else{NULL},
  spring = if(length(spring_out) > 0){do.call("rbind", spring_out)}else{NULL},
  summ = if(length(summ_out) > 0){do.call("rbind", summ_out)}else{NULL},
  fall = if(length(fall_out) > 0){do.call("rbind", fall_out)}else{NULL}
)

#add sp name to the df, catch any null dfs
sp_ls[[j]] <- do.call("rbind", ind_ls) %>%
  {if(!is.null(.)){
    mutate(., species = species$taxon_canonical_name[s])
  }else{.}}

} #j (end loop through individuals)