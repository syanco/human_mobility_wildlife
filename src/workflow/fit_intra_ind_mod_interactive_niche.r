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
  
  .cores <- 4
  .iter <- 1000
  .thin <- 1
  
} else {
  library(docopt)
  # library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'src/funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  # .nc <- ag$nc
  
  # ag <- docopt(doc, version = '0.1\n')
  
  .wd <- getwd()
  .cores <- ag$cores
  .iter  <- ag$iter
  .thin <- ag$thin
  
  source(file.path(.wd,'src/funs/input_parse.r'))
  
  .datPF <- makePath(ag$dat)
  .outP <- makePath(ag$out)
  
}

#---- Initialize Environment ----#
# Start time
t0 <- Sys.time()

# Run startup
source(file.path(.wd,'src/startup.r'))

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
conflict_prefer("has_name", "tibble")

.cores <- ifelse(is.null(.cores), 10, as.numeric(.cores))
.iter <- ifelse(is.null(.iter), 3000, as.numeric(.iter))
.thin <- ifelse(is.null(.thin), 4, as.numeric(.thin))

#Source all files in the auto load funs directory
list.files(file.path(.wd,'src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Load Data ----#

message("Loading data...")

# Load trsait data
traits <- read_csv("raw_data/covid_movement_full_repo/raw_data/anthropause_data_sheet.csv")

# load size data
size <- read_csv("out/dbbmm_size.csv") %>%
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>% 
  mutate(ind_f = as.factor(ind_id))%>%  # create factor version of ind for REs)
  mutate(species = case_when( # correct species names
    study_id == 1442516400 ~ "Anser caerulescens",
    study_id == 1233029719 ~ "Odocoileus virginianus",
    study_id == 1631574074 ~ "Ursus americanus",
    study_id == 1418296656 ~ "Numenius americanus",
    study_id == 474651680  ~ "Odocoileus virginianus",
    study_id == 1044238185 ~ "Alces alces",
    TRUE ~ species
  ))%>% 
  mutate(species = case_when(
    species == "Chen caerulescens" ~ "Anser caerulescens",
    TRUE ~ species
  )) %>% 
  distinct()

# load and process data
breadth <- read_csv("out/niche_determinant_anthropause.csv") %>%
  mutate(scientificname = case_when( # correct species names
    studyid == 1442516400 ~ "Anser caerulescens",
    studyid == 1233029719 ~ "Odocoileus virginianus",
    studyid == 1631574074 ~ "Ursus americanus",
    studyid == 1418296656 ~ "Numenius americanus",
    studyid == 474651680  ~ "Odocoileus virginianus",
    studyid == 1044238185 ~ "Alces alces",
    TRUE ~ scientificname
  ))%>% 
  mutate(scientificname = case_when(
    scientificname == "Chen caerulescens" ~ "Anser caerulescens",
    TRUE ~ scientificname
  )) %>% 
  distinct() %>%
  mutate(tmax_mvnh = tmax,
         # tmin_mvnh = tmin,
         # elev_mvnh = elev,
         ndvi_mvnh = ndvi) %>%
  select(!c(tmax, ndvi)) %>%
  filter(studyid != 351564596) %>%
  filter(studyid != 1891587670) %>%
  mutate(ind_f = as.factor(individual)) %>%  # create factor version of ind for REs)
  left_join(size, by = c("scientificname" = "species", 
                           # "individual" = "ind_id", 
                         "ind_f" = "ind_f",
                           "year" = "year", 
                           "studyid" = "study_id", 
                           "week" = "wk"))


message("Processing the data to allow the magic to happen...")

# identify paired observations only
paired <- breadth %>% 
  # mutate(yrnum = as.numeric(yr)) %>% 
  group_by(individual) %>% 
  summarize(minyr = min(year),
            maxyr = max(year),
            paired = minyr != maxyr) %>% 
  filter(paired)

# create vector of paired individual ids
paired_vec <- paired %>% 
  pull(individual)

# filter size dataset to only inlcude paired individuals
breadth_paired <- breadth %>% 
  filter(individual %in% paired_vec) %>% 
  #filter to mammals only
  left_join(traits, by = c("scientificname" = "Species")) %>% 
  # filter(Family == "Cervidae" | Family == "Canidae" | Family == "Ursidae" |
  #          Family == "Felidae" | Family == "Antilocapridae" | Family == "Bovidae") %>% 
  mutate(diet = case_when(Diet.PlantO >= 50 ~ "Herbivore",
                          Diet.Fruit >= 50 ~ "Frugivore",
                          Diet.Scav >= 50 ~ "Savenger",
                          Diet.Inv >= 50 ~ "Insectivore",
                          (Diet.Vend+Diet.Vect+Diet.Vfish) >= 50 ~ "Carnivore",
                          TRUE ~ "Omnivore"),
         breadth_scale = scale(log(total+0.00000000000001)),
         log_area = log(area), #get log of weekly area use
         log_area_scale = scale(log_area), # standardize it
         sg_norm = sg / cbg_area, # normalize safegraph data by size of the CBG
         sg_sqrt = sqrt(sg_norm),
         sg_scale = scale(sg_norm),
         # log_sg_norm = log(sg_norm),
         ghm_scale = scale(ghm),
         ndvi_scale = scale(ndvi),
         # lst_scale = scale(lst),
         tmax_scale = scale(tmax),
         ind_f = as.factor(individual), # create factor version of ind for REs
         grp = paste(ind_f, year, sep = "_"), # create indXyr grouping factor
         # trt_new = gsub('_.*','',trt),
         year_f = factor(year), # create year factor
         # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
         # sp2 = gsub(" ", "_", species),
         wk_n = as.numeric(substring(week, 2)), # extract week number
         ts = parse_date_time(paste(year, week, 01, sep = "-"), "%Y-%U-%u"), # make better date format
         study_f = as.factor(studyid),
         ind_wk = paste0(individual,week))  # make study id factor) %>% 

# create wide format and calculate deltas
breadth_wide <- breadth_paired %>% 
  pivot_wider(id_cols = c(ind_f, week, scientificname), 
              values_from = c(breadth_scale, log_area_scale, sg_norm, ghm_scale, ndvi_scale, tmax_scale), 
              names_from = year_f) %>% 
  mutate(breadth_diff = breadth_scale_2019-breadth_scale_2020,
         area_diff = log_area_scale_2019-log_area_scale_2020,
         sg_diff = sg_norm_2019-sg_norm_2020,
         ghm_diff = ghm_scale_2019-ghm_scale_2020,
         ndvi_diff = ndvi_scale_2019-ndvi_scale_2020,
         tmax_diff = tmax_scale_2019-tmax_scale_2020) %>% 
  filter(!is.nan(sg_diff)) %>% 
  filter(!is.na(sg_diff))

# # get sample size per sp
# size_wide %>% 
#   group_by(species) %>% 
#   summarize(nind = n_distinct(ind_id)) %>% 
#   arrange(-nind)



#---- Load Data ----#

message("Starting that modeling magic...")

form <-  bf(breadth_diff ~ area_diff + sg_diff*ghm_diff + ndvi_diff + tmax_diff + (1|scientificname/ind_f) + ar(time = week, gr = ind_f))
message("Fitting models with formula:")
print(form)

message("Starting model...")


# fit model
mod <- brm(
  form,
  data = breadth_wide,
  family = student(),
  inits = 0,
  cores = 4,
  iter = 10000,
  thin = 5
)

#stash results into named list
out <- list(
  data = breadth_wide,
  model = mod
)

#write out results
save(out, file = glue("{.outP}/niche_intra_ind_int_mod_{Sys.Date()}.rdata"))


#---- Finalize script ----#


message(glue('Script complete in {diffmin(t0)} minutes'))
