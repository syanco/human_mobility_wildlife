# wd=~/projects/Anthropause/analysis
# src=~/projects/Anthropause/src
# db=data/mosey_mod.db
# 
# 
# cd $wd
# 
# #--- Add a column to the event table to hold the env data
# 
# sql="alter table event add column population_d real"
# 
# sqlite3 $db "$sql;"
# 
# sql="alter table event add column stringency real"
# 
# sqlite3 $db "$sql;"
# 
# sql="alter table event add column county_fips real"
# 
# sqlite3 $db "$sql;"
# 
# sql="alter table event add column census_tract real"
# 
# sqlite3 $db "$sql;"

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 
# Annotate Human population density for ndividual animals
# Next script HFI
# Next script google mobility


# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

#!/usr/bin/env Rscript --vanilla

# This script implements the breezy philosophy: github.com/benscarlson/breezy

# ==== Breezy setup ====

#TODO: Update to allow annotation of multiple layers?

'
Updates the event table with environmental annotations

Usage:
poc-anno_local.r [--db=<db>] [--entity=<entity>] [--seed=<seed>] [-b] [-t]
poc-anno_local.r (-h | --help)

Control files:
ctfs/study.csv
ctfs/individual.csv (optional)


Options:
-h --help     Show this screen.
-v --version     Show version.
-d --db=<db> Path to movement database. Defaults to <wd>/data/move.db
-e --entity=<entity>  Either <study|individual> Defines control file that will be used to run ctmm. Defaults to study
-s --seed=<seed>  Random seed. Defaults to 5326 if not passed
-b --rollback
-t --test         Indicates script is a test run, will not save output parameters or commit to git
' -> doc

#---- Input Parameters ----#

if(interactive()){
  library(here)
  
  .wd <- '~/projects/Anthropause/analysis'
  .seed <- NULL
  .rollback <- TRUE
  .test <- TRUE
  rd <- here::here
  
  # Path to Moset DB
  .dbPF <- '~/projects/Anthropause/analysis/data/mosey_mod.db'
  .entity <- 'study' # setting to study 
  #  .layerPF <- '~/projects/Anthropause/analysis/data/env_layers/HFI_wgs_2013_USA.tif'
  
  # Load oxford stringency:
  require(lubridate)
  oxford_str <- read.csv('/Users/diegoellis/Downloads/USA-covid-policy-master/data/OxCGRT_US_subnational_05Aug2020.csv') %>% select(CountryCode, RegionName, Date, StringencyIndex) %>%
    mutate(date = ymd(Date), state = RegionName) %>% select(-CountryCode, -RegionName, -Date)
  # oxford_str$date <- ymd(oxford_str$Date) %>% mutate(state = RegionName) %>% select(-CountryCode, RegionName)
  # oxford_str$state <- oxford_str$RegionName
  # Load google mobility data:
  # US_google_mobility <- read.csv('/Users/diegoellis/projects/COVID/COVID/Human_mobility_data/Region_Mobility_Report_CSVs/2020_US_Region_Mobility_Report.csv')  %>% select(-iso_3166_2_code, metro_area) %>% drop_na(sub_region_1, sub_region_2)
  # # US_google_mobility <- read.csv('/Users/diegoellis/Downloads/Region_Mobility_Report_CSVs/2020_US_Region_Mobility_Report.csv') %>% select(-iso_3166_2_code, metro_area) %>% drop_na(sub_region_1, sub_region_2)
  # US_google_mobility$date <- ymd(US_google_mobility$date)
  # # sub_region_2 and sub_region 1
  # US_google_mobility[US_google_mobility$sub_region_1=="",]<-NA
  # US_google_mobility[US_google_mobility$sub_region_2=="",]<-NA
  # US_google_mobility = US_google_mobility %>% drop_na(sub_region_1, sub_region_2)
  # US_google_mobility$NAME <- paste0(US_google_mobility$sub_region_2,', ', US_google_mobility$sub_region_1)
  # US_google_mobility$USA_county = gsub(",.*","",US_google_mobility$NAME)
  
  # G_mob = US_google_mobility %>% mutate(NAME = NAME,
  #                                       state = sub_region_1,
  #                                       retail_recreation = retail_and_recreation_percent_change_from_baseline,
  #                                       grocery_pharma = grocery_and_pharmacy_percent_change_from_baseline,
  #                                       greenspace = parks_percent_change_from_baseline,
  #                                       transit = transit_stations_percent_change_from_baseline,
  #                                       workplace = workplaces_percent_change_from_baseline,
  #                                       residential = residential_percent_change_from_baseline
  # ) %>% select(NAME, state, date, retail_recreation, grocery_pharma, greenspace, transit, workplace, residential)
  
  
  # Remove everything after the coma:
  # G_mob$county = gsub(",.*","",G_mob$NAME)
  
  # Select the same time line:
  # lockdown_start <- as.POSIXct('2020-03-01 00:00:00')
  # lockdown_end <- as.POSIXct('2020-06-01 00:00:00')
  # baseline_start <- as.POSIXct('2019-03-10 00:00:00')
  # baseland_end <- as.POSIXct('2019-06-20 00:00:00')
  # 
  # US_google_mobility <- droplevels(subset(US_google_mobility, date >=lockdown_start & date <= lockdown_end))
  # oxford_str <- droplevels(subset(oxford_str, date >=lockdown_start & date <= lockdown_end))
  
  
  # us_stingency_google <- left_join(US_google_mobility, oxford_str)
  # us_states <- left_join(us_states, us_stingency_google)
  
  
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
  .dbPF <- makePath(ifelse(length(ag$db)!=0, ag$db,'data/mosey_mod.db'))
  .entity <- ifelse(is.null(ag$entity),'entity',ag$entity)
  # .layerPF <- ag$layer
}





#---- Initialize Environment ----#
if(!is.null(.seed)) {message(paste('Random seed set to',.seed)); set.seed(as.numeric(.seed))}

t0 <- Sys.time()

source(rd('src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    # Put all your libraries in ehre: 
    library(DBI)
    library(RSQLite)
    library(terra)
    library("rlang", lib.loc=.libPaths())
    require(tidycensus)
    options(tigris_use_cache = TRUE) # speeds future usage up
    require(spData)
    require(sf)
    library(rgeos)
    require(move)
    require(lubridate)
    data(us_states)
    data("fips_codes")
    fips_codes = fips_codes %>% mutate(state = state_code) %>% select(state, state_name)
    
    library(sp)
    library(maps)
    library(maptools)
    conflict_prefer("map", "purrr")
    # The single argument to this function, pointsDF, is a data.frame in which:
    #   - column 1 contains the longitude in degrees (negative in the US)
    #   - column 2 contains the latitude in degrees
    
    latlong2county <- function(pointsDF) {
      # Prepare SpatialPolygons object with one SpatialPolygon
      # per county
      counties <- map('county', fill=TRUE, col="transparent", plot=FALSE)
      IDs <- sapply(strsplit(counties$names, ":"), function(x) x[1])
      counties_sp <- map2SpatialPolygons(counties, IDs=IDs,
                                         proj4string=CRS("+proj=longlat +datum=WGS84"))
      
      # Convert pointsDF to a SpatialPoints object 
      pointsSP <- SpatialPoints(pointsDF, 
                                proj4string=CRS("+proj=longlat +datum=WGS84"))
      
      # Use 'over' to get _indices_ of the Polygons object containing each point 
      indices <- over(pointsSP, counties_sp)
      
      # Return the county names of the Polygons object containing each point
      countyNames <- sapply(counties_sp@polygons, function(x) x@ID)
      countyNames[indices]
    }
    
    
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



#---- Perform analysis ----#

# Start with 1 study 
# i <- 1

for(i in 1:nrow(studies)) { # i <- 1
  #i <- 1
  study <- studies[i,]
  
  message('**********')
  message(glue('Starting annotation for study: {study$study_name} (id: {study$study_id})'))
  

  if(.entity=='individual') {
    inddat <- inds %>% filter(study_id==study$study_id)
  } else {
    inddat <- 'select individual_id, local_identifier 
    from individual where study_id = {study$study_id}' %>%
      glue_sql(.con=db) %>%
      dbGetQuery(db,.) %>% as_tibble
  }
  
  # For every study, for every individuals grab its bounding box, check what USA states are in there and get census tract level informaiton for this
  
  #j <- 1
  for(j in 1:nrow(inddat)) {
    
    
    
    ind <- inddat[j,]
    
    message(glue('Starting annotation for {ind$local_identifier} (id: {ind$individual_id})'))
    
    sql <- 'select event_id,lon,lat,timestamp
    from event
    where individual_id = {ind$individual_id}' %>% glue_sql(.con=db)
    
    message('Getting points from the database')
    tic()
    dat <- dbGetQuery(db, sql) %>% as_tibble() %>% st_as_sf(coords = c("lon", "lat"), crs="+proj=longlat +datum=WGS84")
    toc()
    message(glue('Got {nrow(dat)} points'))
    
    # pts <- vect(dat,crs="+proj=longlat +datum=WGS84") # Terra object 
    pts <- dat 
    
    message('Performing annotation')
    
    # Search for variables of interest. ACS: https://data.census.gov/cedsci/
    # Convert individual to a sf object
    pts = pts %>% st_transform(st_crs(us_states))
    # Intersect the animal bounding box with all US states
    states_of_interest <- st_crop(us_states, st_bbox(pts))
    states_of_interest$state_abbrev <- usdata::state2abbr(states_of_interest$NAME)
    
    
    pts = st_intersection(pts, states_of_interest['NAME']) %>% mutate(state = NAME) %>% select(-NAME)
    
    
    # pts_tmp <- as(as(pts, 'Spatial'), 'data.frame')
    
    # Left join oxford government response trakcer:
    pts$date <- ymd_hms(pts$timestamp)
    pts$date <- as.Date(pts$date)
    
    pts = left_join(pts, oxford_str, by = c('state', 'date'))
    
    
    # Load google stringency index:
    
    # Get the abbreviation of the states in which our animal occurs:
    # Get decenial for census block level ! ###### # Need ot update the variable here ! 
    pop_density <- get_acs(geography = "tract",
                           variables = c(pop_count = "B01003_001"),
                           state = states_of_interest$state_abbrev,
                           cache = TRUE,
                           geometry = TRUE) %>% mutate(human_population = estimate)
    pop_density$Area <- st_area(pop_density)
    
    pop_density$p_density = pop_density$estimate / pop_density$Area
    
    counties_of_interest <- st_crop(pop_density, st_bbox(pts))
    # counties_of_interest <- st_intersection(pop_density, pts)
    
    
    counties_of_interest$NAME = sub(".*?,", "", counties_of_interest$NAME)
    
    counties_of_interest$state <- substr(start = 1, stop = 2, counties_of_interest$GEOID)
    counties_of_interest$county_fips <- substr(start = 1, stop = 5, counties_of_interest$GEOID)
    counties_of_interest$census_tract <- substr(start = 1, stop = 11, counties_of_interest$GEOID)
    # Remove everything after the coma:
    counties_of_interest$NAME = gsub(",.*","",counties_of_interest$NAME)
    counties_of_interest$NAME = substring(counties_of_interest$NAME, 2)
    
    
    counties_of_interest =  left_join(counties_of_interest, fips_codes, by = 'state') %>% select(-state) %>% mutate(state = state_name) %>% select(-state_name)
    
    
    # counties_of_interest$NAME = gsub(".*,","",counties_of_interest$NAME)
    # regmatches(counties_of_interest$NAME,gregexpr("(?<=,).*",counties_of_interest$NAME,perl=TRUE))
    
    # Remove eferything after the comma:
    # G_mob_tmp = G_mob[G_mob$state %in% counties_of_interest$state,]
    # G_mob_tmp = G_mob[G_mob$county %in% counties_of_interest$NAME,]
    # US_google_mobility$USA_county = gsub(",.*","",US_google_mobility$NAME)
    # US_google_mobility[US_google_mobility$USA_county %in% counties_of_interest$NAME, ]
    
    
    # pts = left_join(pts, G_mob, by = c('state', 'date'))
    
    # Annotate static variables of interest:
    pts = pts %>% mutate(Stringency = StringencyIndex) %>% select(-StringencyIndex)
    
    #    ano = st_join(pts, pop_density)
    # ano = st_join(pts, pop_density) %>% st_transform(crs('+proj=longlat +datum=WGS84'))  %>% as_tibble() %>% select(event_id, p_density, Stringency) # if i dont do roads uncomment this 
    
    # a = st_intersection(pts, counties_of_interest)
    # pts = st_intersection(pts, counties_of_interest) %>% mutate(state = NAME) %>% select(-NAME)
    
    
    
    ano = st_join(pts, counties_of_interest) %>% st_transform(crs('+proj=longlat +datum=WGS84')) %>% as_tibble() 
    ano = ano %>% mutate(state = state.x) %>% select(event_id, GEOID, p_density, Stringency, state, county_fips, census_tract)   # if i dont do roads uncomment this 
    ano = ano %>% select(-state, -GEOID)
    
    
    ano = plyr::ddply(ano, 'event_id', function(x){
      data.frame(
      p_density = unique(x$p_density),
      Stringency = unique(x$Stringency),
      county_fips = unique(x$county_fips),
      census_tract = unique(x$census_tract)
      ) %>% as_tibble()
    })
    
    # ano %>% group_by(event_id) %>% summarise()
    
    # ADD GOOGLE MOBILITY DATA
    #     left_join(ano, counties_of_interest, by = c('state', ))
    
    #Initial tests show might be faster if cropping first. Need to experiment more
    #Seems to be true for small extents, but crop is slow for large extents
    #rst <- rast(env$path) %>% crop(ext(pts))
    # rst <- rast(.layerPF)
    
    # ano <- terra::extract(rst,pts) %>% as_tibble #1e5, 0.283 sec
    toc()
    
    #Need to rely on correct row order to join back points. Feels dangerous!
    # Also, can't use name, so use position. Should always be ID,<layer name>
    # udat <- pts %>% 
    #   as.data.frame %>% as_tibble %>%
    #   mutate(!!.colname := ano %>% pull(2)) %>%
    #   select(event_id,!!.colname)
    # 
    #-- Update the database
    sql <- 'update event set population_d = $p_density, stringency=$Stringency, county_fips=$county_fips, census_tract=$census_tract
     where event_id = $event_id' %>% glue
    
    
    
    # # -- Update the database
    # sql <- 'update event set Stringency = $Stringency, Park_m_G = $parks,
    #  Workplace_m_G = $workplace, Recreation_m_G = $recreation_mob, 
    #  Grocery_m_G = $grocery, Residential_m_G = $residential
    #   where event_id = $event_id' %>% glue
    # 
    #  sql <- 'update event set population_d = $p_density
    # where event_id = $event_id' %>% glue
    
    dbExecute(db,'PRAGMA foreign_keys=ON')
    dbBegin(db) # Here the transaction begins 
    
    message(glue('Updating event table'))
    tic()
    rs <- dbSendStatement(db, sql) #parameter names should match column names
    # Ask Ben for help here #####
    dbBind(rs,params=ano) #just pass in the full dataframe #ANO is my data
    rows <- dbGetRowsAffected(rs)
    dbClearResult(rs)
    toc()
    
    #Close transaction
    if(.rollback) {
      message('Rolling back transaction because this is a test run.')
      dbRollback(db)
    } else {
      if(rows==nrow(dat)) {dbCommit(db)} else {message('Rolling back'); dbRollback(db)} # db rollback, undo everything oyu made. 
    }
  }
}
#---- Finalize script ----#

if(!.test) {
  #  library(git2r)
  # library(uuid)
  
  .runid <- UUIDgenerate()
  .parPF <- file.path(.wd,"run_params.csv")
  
  #Update repo and pull out commit sha
  repo <- repository(rd('src'))
  
  rstat <- status(repo)
  if(length(rstat$staged) + 
     length(rstat$unstaged) + 
     length(rstat$untracked) > 0) {
    add(repo,'.')
    commit(repo, glue('script auto update. runid: {.runid}'))
  }
  
  
  .git_sha <- sha(repository_head(repo))
  
  #Save all parameters to csv for reproducibility
  #TODO: write this to a workflow database instead
  saveParams(.parPF)
}

dbDisconnect(db)

message(glue('Script complete in {diffmin(t0)} minutes'))
