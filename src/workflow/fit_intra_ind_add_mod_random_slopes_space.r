#!/usr/bin/env Rscript 
#
# This script models random slopes within species as a function of changes in human activities
# and environmental conditions experienced by each individual with movement data in both 2019 
# and 2020. Human mobility and landscape modification are considered additive terms. The 
# model includes NDVI and TMAX as additive fixed effects, random slopes by species, and groups 
# the individual random intercept hierarchically within species. Afterwards, we check
# autocorrelation with the script check_intra_ind_mod_ac.R
 
# See manuscript section: Behavioral plasticity of responses

# ==== Setup ====

'
Usage:
fit_intra_ind_mod_random_slopes.r <dat> <out> <cores> [<iter> <thin>]
fit_intra_ind_mod_random_slopes.r (-h | --help)

Parameters:
  dat: path to input data
  out: path to output directory.
  cores: number of cores for parallel processing
  
Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  
  .wd <- '/home/julietcohen/repositories/human_mobility_wildlife'
  
  .datPF <- file.path(.wd,'out/dbbmm_size.csv')
  .outP <- file.path(.wd,'out/intra_ind_models')
  
  .cores <- 1
  .iter <- 10000
  .thin <- 5
  
} else {
  library(docopt)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'src/funs/input_parse.r'))
  
  .datPF <- makePath(ag$dat)
  .outP <- makePath(ag$out)
  .cores <- ag$cores
  .iter  <- ag$iter
  .thin <- ag$thin
  
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

# Load trait data
traits <- read_csv("/home/julietcohen/covid_movement_full_repo/raw_data/anthropause_data_sheet.csv") %>%
  add_row(Species = "Procyon lotor",
          class = "mammal",
          migratory = "non-migratory") %>% 
  add_row(Species = "Spilogale putorius",
          class = "mammal",
          migratory = "non-migratory") %>% 
  add_row(Species = "Sus scrofa",
          class = "mammal",
          migratory = "non-migratory") %>% 
  mutate(Species = case_when(Species == "Chen rossii" ~ "Anser rossii",
                             TRUE ~ Species))

# load size data
size <- read_csv(file.path(.datPF)) %>%
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>% 
  mutate(ind_f = as.factor(ind_id))%>%  # create factor version of ind for REs)
  mutate(species = case_when( # correct species names
    study_id == 1442516400 ~ "Anser caerulescens",
    study_id == 1631574074 ~ "Ursus americanus",
    study_id == 1418296656 ~ "Numenius americanus",
    study_id == 474651680  ~ "Odocoileus virginianus",
    study_id == 1044238185 ~ "Alces alces",
    TRUE ~ species
  ))%>% 
  mutate(species = case_when(
    species == "Chen caerulescens" ~ "Anser caerulescens",
    species == "Chen rossii" ~ "Anser rossii",
    TRUE ~ species
  )) %>% 
  distinct()


message("Processing the data to allow the magic to happen...")

# identify paired observations only
paired <- size %>% 
  # mutate(yrnum = as.numeric(yr)) %>% 
  group_by(ind_f) %>% 
  summarize(minyr = min(year),
            maxyr = max(year),
            paired = minyr != maxyr) %>% 
  filter(paired) %>% 
  pull(ind_f)

# filter size dataset to only inlcude paired individuals
size_paired <- size %>% 
  filter(ind_f %in% paired) %>% 
  #filter to mammals only
  left_join(traits, by = c("species" = "Species")) %>% 
  # filter(Family == "Cervidae" | Family == "Canidae" | Family == "Ursidae" |
  #          Family == "Felidae" | Family == "Antilocapridae" | Family == "Bovidae") %>% 
  mutate(diet = case_when(Diet.PlantO >= 50 ~ "Herbivore",
                          Diet.Fruit >= 50 ~ "Frugivore",
                          Diet.Scav >= 50 ~ "Savenger",
                          Diet.Inv >= 50 ~ "Insectivore",
                          (Diet.Vend+Diet.Vect+Diet.Vfish) >= 50 ~ "Carnivore",
                          TRUE ~ "Omnivore"),
         # breadth_scale = scale(log(total+0.000000000000000000001)),
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
         # ind_f = as.factor(individual), # create factor version of ind for REs
         grp = paste(ind_f, year, sep = "_"), # create indXyr grouping factor
         # trt_new = gsub('_.*','',trt),
         year_f = factor(year), # create year factor
         # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
         # sp2 = gsub(" ", "_", species),
         # wk_n = as.numeric(substring(week, 2)), # extract week number
         # ts = parse_date_time(paste(year, week, 01, sep = "-"), "%Y-%U-%u"), # make better date format
         study_f = as.factor(study_id),
         ind_wk = paste0(ind_f,wk))  # make study id factor) %>% 

# create wide format and calculate deltas
size_wide <- size_paired %>% 
  pivot_wider(id_cols = c(ind_f, wk, species), 
              values_from = c(log_area_scale, sg_norm, ghm_scale, ndvi_scale, tmax_scale), 
              names_from = year_f) %>% 
  mutate(size_diff = log_area_scale_2019-log_area_scale_2020,
         sg_diff = sg_norm_2019-sg_norm_2020,
         ghm_diff = ghm_scale_2019-ghm_scale_2020,
         ndvi_diff = ndvi_scale_2019-ndvi_scale_2020,
         tmax_diff = tmax_scale_2019-tmax_scale_2020) %>% 
  filter(!is.nan(sg_diff)) %>% 
  filter(!is.na(sg_diff)) %>%
  # sort data by week within ind within sp 
  arrange(species, ind_f, wk)

# count number of individuals per species in the paired data
# filter to 5+
spp_sufficient_ss <- size_wide %>% 
                    group_by(species) %>% 
                    summarise(n_ind = n_distinct(ind_f)) %>%
                    filter(n_ind >= 5) 

# subset paired data to just the species with 5+ individuals
size_wide_sub <- size_wide %>% 
  semi_join(spp_sufficient_ss, by = "species") %>%
  # scale each difference value 
  ungroup() %>%
  mutate(size_diff = scale(size_diff, center = F),
         ghm_diff = scale(ghm_diff, center = F),
         sg_diff = scale(sg_diff, center = F),
         ndvi_diff = scale(ndvi_diff, center = F),
         tmax_diff = scale(tmax_diff, center = F)) %>%
  # filter outliers for safegraph values
  filter(sg_diff < 2.5 & sg_diff > -2.5)


#---- Load Data ----#

message("Starting that modeling magic...")
# form <- bf(size_diff ~ 1 + sg_diff + ghm_diff + ndvi_diff + tmax_diff + (1 + sg_diff + ghm_diff + ndvi_diff + tmax_diff | species) + (1 | species:ind_f)) 
form <- bf(size_diff ~ 1 + sg_diff + ghm_diff + ndvi_diff + tmax_diff + (1 + sg_diff + ghm_diff + ndvi_diff + tmax_diff || species) + (1 | gr(ind_f, by = species)))
message("Fitting models with formula:")
print(form)

message("Starting model...")


# fit model
mod <- brm(
  form,
  data = size_wide_sub,
  family = student(),
  init = 0,
  cores = .cores,
  iter = .iter,
  thin = .thin,
  control = list(adapt_delta = 0.95)
)

#stash results into named list
out <- list(
  data = size_wide_sub,
  model = mod
)

#write out results
save(out, file = glue("{.outP}/size_intra_ind_add_rs_mod_{Sys.Date()}.rdata"))

#---- Finalize script ----#

message(glue('Script complete in {diffmin(t0)} minutes'))
