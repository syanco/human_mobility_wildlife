#!/usr/bin/env Rscript 
#
#CONDA: covid
#
# This script generates individual dynamic brownian bridge models and associated 
# UDs for migratory  periods.

# TODO:  The dBBMM paramaters (e.g., window size, margin, error, etc.) are 
# currently hardcoded.  Could be passed in as options to the script.
# TODO: verify the volume/probability problem for write out.
# 
# ==== Setup ====

'
Usage:
fit_intra_ind_mod_additive.r <dat> <out> <cores> [<iter> <thin>]
fit_intra_ind_mod_additive.r (-h | --help)

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
  library(here)
  
  .wd <- getwd()
  
  .datPF <- file.path(.wd,'out/dbbmm_size.csv')
  .outP <- file.path(.wd,'out/intra_ind_models')
  
  .cores <- 20
  .iter <- 5000
  .thin <- 4
  
} else {
  library(docopt)
  # library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .wd <- getwd()
  .cores <- ag$cores
  .iter  <- ag$iter
  .thin <- ag$thin
  
  .datPF <- makePath(ag$dat)
  .outP <- makePath(ag$out)
  
}

#---- Initialize Environment ----#
# Start time
t0 <- Sys.time()

# check arg inputs
.cores <- ifelse(is.null(.cores), 4, as.numeric(.cores))
.iter <- ifelse(is.null(.iter), 10000, as.numeric(.iter))
.thin <- ifelse(is.null(.thin), 5, as.numeric(.thin))

# Run startup
source(file.path(.wd,'analysis/src/startup.r'))

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    # library(DBI)
    # library(RSQLite)
    library(lubridate)
    # library(raster)
    # library(move)
    library(doMC)
    library(foreach)
    library(brms)
  }))

conflict_prefer("accumulate", "purrr")
conflict_prefer("ar", "brms")
conflict_prefer("lag", "stats")
conflict_prefer("when", "purrr")

#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Load Data ----#

message("Loading data...")

# Load trsait data
traits <- read_csv("raw_data/anthropause_data_sheet.csv")

# load size data
size <- read_csv("out/dbbmm_size.csv") %>%
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>%
  distinct()

message("Processing the data to allow the magic to happen...")

# identify paired observations only
paired <- size %>% 
  # mutate(yrnum = as.numeric(yr)) %>% 
  group_by(ind_id) %>% 
  summarize(minyr = min(year),
            maxyr = max(year),
            paired = minyr != maxyr) %>% 
  filter(paired)

# create vector of paired individual ids
paired_vec <- paired %>% 
  pull(ind_id)

# filter size dataset to only inlcude paired individuals
size_paired <- size %>% 
  filter(ind_id %in% paired_vec) %>% 
  #filter to mammals only
  left_join(traits, by = c("species" = "Species")) %>% 
  filter(Family == "Cervidae" | Family == "Canidae" | Family == "Ursidae" |
           Family == "Felidae" | Family == "Antilocapridae" | Family == "Bovidae") %>% 
  mutate(diet = case_when(Diet.PlantO >= 50 ~ "Herbivore",
                          Diet.Fruit >= 50 ~ "Frugivore",
                          Diet.Scav >= 50 ~ "Savenger",
                          Diet.Inv >= 50 ~ "Insectivore",
                          (Diet.Vend+Diet.Vect+Diet.Vfish) >= 50 ~ "Carnivore",
                          TRUE ~ "Omnivore"),
         log_area = log(area), #get log of weekly area use
         log_area_scale = scale(log_area), # standardize it
         sg_norm = sg / cbg_area, # normalize safegraph data by size of the CBG
         sg_sqrt = sqrt(sg_norm),
         sg_scale = scale(sg_norm),
         # log_sg_norm = log(sg_norm),
         ghm_scale = scale(ghm),
         ndvi_scale = scale(ndvi),
         lst_scale = scale(lst),
         ind_f = as.factor(ind_id), # create factor version of ind for REs
         grp = paste(ind_f, year, sep = "_"), # create indXyr grouping factor
         # trt_new = gsub('_.*','',trt),
         year_f = factor(year), # create year factor
         # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
         # sp2 = gsub(" ", "_", species),
         wk_n = as.numeric(substring(wk, 2)), # extract week number
         ts = parse_date_time(paste(year, wk, 01, sep = "-"), "%Y-%U-%u"), # make better date format
         study_f = as.factor(study_id),
         ind_wk = paste0(ind_id,wk))  # make study id factor) %>% 

# create wide format and calculate deltas
size_wide <- size_paired %>% 
  pivot_wider(id_cols = c(ind_id, wk, species), 
              values_from = c(log_area_scale, sg_norm, ghm_scale), 
              names_from = year_f) %>% 
  mutate(area_diff = log_area_scale_2020-log_area_scale_2019,
         sg_diff = sg_norm_2020-sg_norm_2019,
         ghm_diff = ghm_scale_2020-ghm_scale_2019) %>% 
  filter(!is.nan(sg_diff)) %>% 
  filter(!is.na(sg_diff))

# # get sample size per sp
# size_wide %>% 
#   group_by(species) %>% 
#   summarize(nind = n_distinct(ind_id)) %>% 
#   arrange(-nind)



#---- Load Data ----#

message("Starting that modeling magic...")

form <-  bf(area_diff ~ 1 + sg_diff + ghm_diff + (1|species/ind_id) + ar(time = wk, gr = ind_id))
message("Fitting models with formula:")
print(form)

message("Starting model...")

# fit model
mod <- brm(
  form,
  data = size_wide,
  # family = Gamma(link = "log"),
  inits = 0,
  cores =.cores,
  iter = .iter,
  thin = 5
)

#stash results into named list
out <- list(
  data = size_wide,
  model = mod
)

#write out results
save(out, file = glue("{.outP}/intra_ind_add_mod_{Sys.Date()}.rdata"))
