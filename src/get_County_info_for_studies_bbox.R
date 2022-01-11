# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
#
# MPYC Anthropause study
# 
# The aim of this script is to identify the US counties in which we have animal movement data
#
# Based on US county information, we can proceed to download the human mobility data for those counties at the census tract level (next script)
#
# diego.ellissoto@yale.edu
#
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# For this script to work, you will need to have an authentificaiton key to use the US census API tigrid package.
# See 5.3 enable census API key: https://geodacenter.github.io/opioid-environment-toolkit/getACSData-tutorial.html


if(interactive()){
  library(here)
  
  .wd <- '~/projects/Anthropause/analysis'
  .seed <- NULL
  .rollback <- TRUE
  .test <- TRUE
  rd <- here::here
  
  .colname <- 'human_footprint'
  .dbPF <- '~/projects/Anthropause/analysis/data/mosey.db'
  .entity <- 'study' # setting to study 
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  .script <-  whereami::thisfile()
  .seed <- ag$seed
  .test <- as.logical(ag$test)
  .rollback <- as.logical(ag$rollback)
  rd <- is_rstudio_project$make_fix_file(.script) # Root directory 
  
  source(rd('src/funs/input_parse.r')) # Looks for a Rproj file and sets the root. 
  
  .colname <- ag$colname
  .dbPF <- makePath(ifelse(length(ag$db)!=0, ag$db,'data/mosey.db'))
  .entity <- ifelse(is.null(ag$entity),'entity',ag$entity)
  # .layerPF <- ag$layer
}

#---- Initialize Environment ----#
if(!is.null(.seed)) {message(paste('Random seed set to',.seed)); set.seed(as.numeric(.seed))}

t0 <- Sys.time()

source(rd('src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    # Put all your libraries in here: 
    library(DBI)
    library(RSQLite)
    library(terra)
    library("rlang", lib.loc=.libPaths())
    require(tidycensus)
    require(spData)
    require(sf)
    library(tigris)
    library(rgeos)
    require(rgeos)
    data(us_states)
    data(alaska)
    # us_states = st_union( us_states , alaska)
    
    # Load your usa census api key
    census_api_key("db9b3481879b9e79eb8c86608656c3c8a8640bbb", install = TRUE, overwrite=TRUE)
    Sys.getenv("CENSUS_API_KEY")
    options(tigris_use_cache = TRUE)
  }))

#Source all files in the auto load funs directory
list.files(rd('src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Local parameters ----#

#---- Load control files ----#
studies <- read_csv(file.path(.wd,'ctfs/study.csv')) %>%
  filter(as.logical(run)) %>% select(-run)

if(.entity=='individual') {
  #Load inds and constrain by study control file.
  inds <- read_csv(file.path(.wd,'ctfs/individual.csv'),col_types=cols()) %>%
    filter(as.logical(run)) %>% select(-run) %>%
    inner_join(studies %>% select(study_id),by='study_id')
  
  #Now need to constrain studies by the selected individuals
  studies <- studies %>% filter(study_id %in% unique(inds$study_id))
}

#---- Initialize database ----#
invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF)
invisible(assert_that(length(dbListTables(db))>0))


stdtb <- tbl(db,'study')
indtb <- tbl(db,'individual')
evttb <- tbl(db,'event')
alaska = alaska %>% st_transform(st_crs(us_states))

# i <- 50
#---- Perform analysis ----#
for(i in 1:nrow(studies)) { # i <- 23
  #i <- 1
  study <- studies[i,]
  
  message('**********')
  message(glue('Starting annotation for study: {study$study_name} (id: {study$study_id})'))
  
  # Get the bounding box of each study we are looping through
  
  query <- 'select min(lon), max(lon), min(lat), max(lat) from event where individual_id in (select individual_id from individual where study_id = {study$study_id} )'  %>%
    glue_sql(.con=db) %>%
    dbGetQuery(db,.) %>% as_tibble
  
  # Create potential error files
  if( all(is.na(query)) == TRUE){
    print(paste0('Wrong bounding box for ', study))
    write.csv(study, file = paste0('/Users/diegoellis/Desktop/FINESST_2022/Unique_census_tract_covid/County_error/county_id_errror_', study$study_id, '.csv'));
    next
  }
  
  # Transform the bounding box into a spatial object with NAD83 projection
  bbox_of_study = st_as_sfc( st_bbox(c(xmin = query$`min(lon)`, xmax = query$`max(lon)`, ymin = query$`min(lat)`, ymax = query$`max(lat)` ), crs = 4326) )  %>%
    st_transform(st_crs(us_states))
  
  # Montheith ata is UTM 12 
  if(grepl('montheith', study$study_name, fixed = TRUE)){
    print('mntheith strikes again')
    if(study$study_name != 'tmp_muledeer_COVID_montheith'){
      bbox_of_study = as( st_as_sfc( st_bbox(c(xmin = query$`min(lon)`, xmax = query$`max(lon)`, ymin = query$`min(lat)`, ymax = query$`max(lat)` ), crs = NA) ),'Spatial')
      proj4string(bbox_of_study) <- '+proj=utm +zone=12 +ellps=WGS84 +units=m +no_defs'
      bbox_of_study = spTransform(bbox_of_study, crs(us_states))
      bbox_of_study = st_as_sfc(bbox_of_study)
    }
  }
  
  message('Performing annotation')
  
  
  # Add sf object of Alaska and merge with lower 48 sf object. Not doing this within this loop led to problems.
  
  # Intersect the animal bounding box with all US states
  states_of_interest <- st_crop(us_states, bbox_of_study)
  states_of_interest_alaska <- st_crop(alaska, bbox_of_study)
  
  if ( nrow(states_of_interest) == 0){
    message(paste0( study$study_name, ' not located in lower 48'))
  }
  if(! nrow(states_of_interest_alaska) == 0){
    message(paste0( study$study_name, ' has data in Alaska !'))
  }
  
  if( ! nrow(states_of_interest) == 0 & ! nrow(states_of_interest_alaska) == 0){
    message(paste0( study$study_name, ' has data in lower 48 AND Alaska !'))
    states_of_interest = rbind(states_of_interest_alaska, states_of_interest)  
  }
  
  if(nrow(states_of_interest) == 0){
    states_of_interest = states_of_interest_alaska
  }
  
  
  states_of_interest$state_abbrev <- usdata::state2abbr(states_of_interest$NAME)
  
  # Sanity check
  if(nrow(states_of_interest)==0){
    message('Data is not present anywhere in the USA!');next}
  
  # Actual query: Get county level information, for variable population density (does not really matter which, as we just want geometry=TRUE in this call to get us the sf objects. Cache = TRUE will make future cashing of this query call faster)
  pop_density <- get_acs(geography = "county",
                         variables = c(pop_count = "B01003_001"),
                         state = states_of_interest$state_abbrev,
                         cache = TRUE,
                         geometry = TRUE) %>% mutate(human_population = estimate)
  
  
  counties_of_interest <- st_crop(pop_density, bbox_of_study)
  
  counties_of_interest$county <- counties_of_interest$GEOID
  counties_of_interest$states <- substr(counties_of_interest$GEOID, 1, 2)
  counties_of_interest$movebank_study_id <- study$study_id
  
  # Save area of interest: 
  counties_of_interest_df = counties_of_interest %>% as.data.frame() %>% select(-geometry, -moe, -estimate) 
  
  write.csv(counties_of_interest_df, file = paste0('/Users/diegoellis/Desktop/FINESST_2022/Unique_census_tract_covid/', unique(counties_of_interest_df$movebank_study_id), '.csv'))
  
  # Add movebank information
  
}

# dir.create('/Users/diegoellis/Desktop/Unique_census_tract_covid/')


dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))



# There are a couple studies missing;

data_170345563_argos<- getMovebankData(study=170345563, login=loginStored, removeDuplicatedTimestamps=T)

# data_170345563_argos_sf = st_as_sf(data_170345563_argos) %>% st_transform(st_crs(us_states))

data_20550947_argos<- getMovebankData(study=20550947, login=loginStored, removeDuplicatedTimestamps=T)

# data_20550947_argos_sf = st_as_sf(data_20550947_argos) %>% st_transform(st_crs(us_states))

data_42451582_argos<- getMovebankData(study=42451582, login=loginStored, removeDuplicatedTimestamps=T)

# data_42451582_argos_sf = st_as_sf(data_42451582_argos) %>% st_transform(st_crs(us_states))

data_1701962556_vhf<- getMovebankData(study=1701962556, login=loginStored, removeDuplicatedTimestamps=T)

# data_1701962556_vhf_sf = st_as_sf(data_1701962556_vhf) %>% st_transform(st_crs(us_states))



montheith<- getMovebankData(study=1891172051, login=loginStored, removeDuplicatedTimestamps=T)


get_county_id_for_missing_studies <- function(movestack_animals){
  data(us_states)
  data(alaska)
  alaska = alaska %>% st_transform(st_crs(us_states))
  
  sf_object_animals = st_as_sf(movestack_animals) %>% st_transform(st_crs(us_states))
  
  bbox_of_study =   st_as_sfc(st_bbox(sf_object_animals))
  # Intersect the animal bounding box with all US states
  states_of_interest <- st_crop(us_states, bbox_of_study)
  # states_of_interest_tmp <- st_crop(us_states, bbox_of_study)
  states_of_interest_alaska <- st_crop(alaska, bbox_of_study)
  
  
  if ( nrow(states_of_interest) == 0){
    message(paste0( study$study_name, ' not located in lower 48'))
  }
  if(! nrow(states_of_interest_alaska) == 0){
    message(paste0( study$study_name, ' has data in Alaska !'))
  }
  
  if( ! nrow(states_of_interest) == 0 & ! nrow(states_of_interest_alaska) == 0){
    message(paste0( study$study_name, ' has data in lower 48 AND Alaska !'))
    states_of_interest = rbind(states_of_interest_alaska, states_of_interest)  
  }
  
  if(nrow(states_of_interest) == 0){
    states_of_interest = states_of_interest_alaska
  }
  
  
  states_of_interest$state_abbrev <- usdata::state2abbr(states_of_interest$NAME)
  
  if(nrow(states_of_interest)==0){
    message('Data is not present anywhere in the USA!');next}
  
  
  pop_density <- get_acs(geography = "county",
                         variables = c(pop_count = "B01003_001"),
                         state = states_of_interest$state_abbrev,
                         cache = TRUE,
                         geometry = TRUE) %>% mutate(human_population = estimate)
  
  
  counties_of_interest <- st_crop(pop_density, bbox_of_study)
  
  counties_of_interest$county <- counties_of_interest$GEOID
  counties_of_interest$states <- substr(counties_of_interest$GEOID, 1, 2)
  
  counties_of_interest$movebank_study_id <- study$study_id
  # counties_of_interest$movebank_individual_id <- ind$individual_id
  # counties_of_interest$unique_id <- paste0(study$study_id, '_', ind$individual_id)
  
  # Save shapefile 
  
  counties_of_interest_df = counties_of_interest %>% as.data.frame() %>% select(-geometry, -moe, -estimate) 
  
  write.csv(counties_of_interest_df, file = paste0('/Users/diegoellis/Desktop/FINESST_2022/Unique_census_tract_covid/', movestack_animals@study, '.csv'))
  
}

get_county_id_for_missing_studies( data_170345563_argos)
get_county_id_for_missing_studies( data_20550947_argos)

get_county_id_for_missing_studies( data_42451582_argos)
get_county_id_for_missing_studies( data_1701962556_vhf)

# Load all .csv files into a single data frame a <- (rbind.row(lappy(list.files(folder_with_csv), read_csv))) %>% distinct(GEOID)

# 
# 
# # Now look at unique GEO ID:
# 
# a = list.files('/Users/diegoellis/Desktop/FINESST_2022/Unique_census_tract_covid/')
# a = a[!a %in% c('County_error')]
# temp = paste0('/Users/diegoellis/Desktop/FINESST_2022/Unique_census_tract_covid/', a)
# yMergedData <- 
#   do.call(rbind,
#           lapply(temp, read.csv))
# 
# 
# myfiles = lapply(temp, read.delim)
#   
