#!/usr/bin/env Rscript 

# This script implements the breezy philosophy: github.com/benscarlson/breezy

# activate conda env prior to running

#TODO: default to using the band name as the output column name of the env variable
#TODO: add optional parameter to not start task (for debugging)
#TODO: handle exception where point/env data can't be found
#TODO: right now, all datasets need to be in the same folder (datP)
# ==== Breezy setup ====

'
Annotates gee point datasets over multiple entities and variables.

Usage:
anno_gee.r <dat> <out> <ent> <env_id> <colname> <band> 
anno_gee.r (-h | --help)

Control files:

Parameters:
  dat: folder containing feature collections to annotate
  out: path to output directory and file on gcs. do not include file extension, url info or bucket
  ent: the focal study id
  env_id: gee id
  colname: what to call the var
  band: band form the gee image
  
Options:
-h --help     Show this screen.
-v --version     Show version.

' -> doc

#---- Input Parameters ----#
if(interactive()) {

  .pd <- '/Users/juliet/Documents/OliverLab/covid_paper/repositories/human_mobility_wildlife'
  .wd <- file.path(.pd,'')
  .seed <- NULL
  rd <- here::here
  .col_name <- "column"
  .env_id <- "NASA/ORNL/DAYMET_V4"
  .band <- 4
  .std <- "481458"
  
  # Required parameters
  .datP <- 'projects/covid-mvmnt-2024-440720/assets/tracks/481458'
  .outP <- 'annotated/481458_column'
  
  # Optional parameters
  .groups <- NULL #only run for these groups. handy for testing
  .npts <- NULL
  
} else {
  suppressWarnings(
    suppressPackageStartupMessages({
        library(docopt)
        library(rprojroot)
  }))

  ag <- docopt(doc, version = '0.1\n')

  .wd <- getwd()
  .script <-  suppressWarnings(thisfile())
  .seed <- ag$seed
  # rd <- is_rstudio_project$make_fix_file(.script)
  
  source('src/funs/input_parse.r')
  
  # Required parameters

  .datP <- ag$dat
  .outP <- ag$out
  .std <- ag$ent
  .env_id <- ag$env_id
  .col_name <- ag$colname
  .band <- as.numeric(ag$band)
  
  # Optional parameters
  
  # ifelse can't return NULL
  if(is.null(ag$groups)) {.groups <- NULL } else {.groups <- as.integer(parseCSL(ag$groups))}
  if(is.null(ag$npts)) {.npts <- NULL} else {.npts <- as.integer(ag$npts)}
    
}

#---- Initialize Environment ----#

message("Initializing environment...")

if(!is.null(.seed)) {message(paste('Random seed set to',.seed)); set.seed(as.numeric(.seed))}

t0 <- Sys.time()

source('src/startup.r')

# Source all files in the auto load funs directory
list.files('src/funs/auto',full.names=TRUE) %>% walk(source)

reticulate::use_condaenv("hmw_py3-9", required = TRUE)
message("Current Python Path: ", Sys.getenv("PYTHONPATH"))

suppressWarnings(
  suppressPackageStartupMessages({
    library(rgee)
  }))

# Initialize gee
ee_Initialize(user = "jscohen@ucsb.edu",
              quiet = FALSE,
              auth_mode = "notebook",
              auth_quiet = TRUE,
              drive = FALSE,
              gcs = TRUE)
  
#TODO: do a check to make sure rgee initialized correctly

#---- Local parameters ----#
.entity <- 'study'
.colImageId <- 'image_id';
.colMillis <- 'millis'
.colTimestamp <- 'timestamp'
.colGrp <- 'grp'
.bucket <- 'covid-mvmnt-2024'

#---- Control files ----#

entities <- read_csv(file.path(.wd,'ctfs',glue('{.entity}.csv'))) %>% filter(run==1)

#---- Perform analysis ----#
  
message(glue('Processing {.datP}'))

pts <- ee$FeatureCollection(file.path(.datP))

message(glue('Setting up annotation tasks for {.col_name}'))

#Check if the layer is a computed layer. If so load it.
#Otherwise load gee asset
if(file.exists(file.path("~", .env_id))) {
  source(file.path("~", .env_id))
  layer <- getLayer()
  assetType <- getAssetType()
} else {
  assetType <- ee$data$getAsset(.env_id)$type
  
  if(assetType=='IMAGE') {
    layer <- ee$Image(.env_id)$select(list(.band))
  } else if (assetType=='IMAGE_COLLECTION') {
    layer <- ee$ImageCollection(.env_id)
  }  else {
    stop(glue('Invalid asset type: {assetType}'))
  }
}

#For testing
# pts <- ee$FeatureCollection(pts$toList(10))
# pts$aggregate_array(.colGrp)$getInfo()


if(is.null(.groups)) {
  #Groups run from 0...n, so to get number of groups need to add 1
  maxgrp <- pts$aggregate_max(.colGrp)$getInfo()
  groups <- 0:maxgrp
} else {
  groups <- .groups
  message(glue('Running only for the following group numbers: {groups}'))
}

if(length(groups)>1) message(glue('Splitting annotation into {length(groups)} tasks'))

# Note groups start at 0
for(group in groups) {
  
  ptsGrp <- pts$filter(ee$Filter$eq(.colGrp,group))

  # Make sure there are points in the group. Can result in 0 records if .group
  # does not exist in the dataset.
  invisible(assert_that(ptsGrp$size()$getInfo() > 0))
  
  if(!is.null(.npts)) {
    message(glue('Running for a subset of {.npts} points'))
    ptsGrp <- ee$FeatureCollection(ptsGrp$toList(.npts))
  }
  

  if(assetType=='IMAGE') {
    #Note that the band is selected above, when loading the layer
    anno <- layer$reduceRegions(
      reducer = ee$Reducer$median()$setOutputs(list(.col_name)),
      collection = ptsGrp,
      scale = layer$projection()$nominalScale()
    )
    
  } else if(assetType=='IMAGE_COLLECTION') {
    
    #env <- dist2water_month()
    #env <- ee$ImageCollection(.envPF)
    
    ptsGrp <- ptsGrp$map(function(f) {
      mil = ee$Date(f$get(.colTimestamp))$millis()
      f <- f$set(.colMillis,mil)
      return(f)
    })
    
    filter <- ee$Filter$And(
      ee$Filter$lessThanOrEquals('system:time_start', NULL, .colMillis),
      ee$Filter$greaterThan('system:time_end', NULL, .colMillis)
    )
    
    joined <- ee$Join$saveAll('features')$apply(layer, ptsGrp, filter)
    
    anno <- joined$map(function(img) {
      
      img <- ee$Image(img)$select(list(.band))
      
      #View(img$projection()$getInfo())
      #img$projection()$nominalScale()$getInfo()
      
      fc <- ee$FeatureCollection(ee$List(img$get('features')))
      
      vals <- img$reduceRegions(
        #TODO: if I'm just extracting the pixel value, should I use
        # ee$Reducer$first() instead ?
        reducer=ee$Reducer$median()$setOutputs(list(.col_name)),
        scale=img$projection()$nominalScale(),
        collection=fc)
      
      vals <- vals$map(function(f) {
        f$set(.colImageId,img$get('system:index'))
      })
      
      return(vals)
    })$flatten()
    
  } else {
    stop(glue('Invalid asset type: {assetType}'))
  }
  
  anno <- anno$sort('anno_id')
  #View(anno$getInfo()); quit()
  
  fileN <- glue('{.std}_{.col_name}_{group}')
  
  task <- ee$batch$Export$table$toCloudStorage(
    collection=anno,
    description=fileN,
    bucket=.bucket,
    fileNamePrefix=file.path(.outP,fileN),
    fileFormat='csv',
    selectors=list('anno_id', .col_name)
  )
  
  task$start()
  message(glue('Task for group {group} started'))
  message(glue('Results will be saved to gs://{.bucket}/{.outP}/{fileN}.csv'))

} #loop over groups

#---- Finalize script ----#
message(glue('Script complete in {diffmin(t0)} minutes'))
