#!/usr/bin/env Rscript --vanilla

# ==== Input parameters ====

'
Usage:
fit-space-use-models.r <dat> <out> <trait> <cores> [<iter> <thin>]
fit-space-use-models.r (-h | --help)


Parameters:
dat: path to input csv file.
out: path to output directory.
trait: path to species trait csv
cores: number of HPC cores 
iter: MCMC iterations
thin: MCMC thin rate

Options:
-h --help           Show this screen.
-v --version        Show version.

' -> doc

if(interactive()) {
  
  .wd <- getwd()
  
  .datPF <- file.path(.wd,'out/dbbmm_size.csv')
  .outP <- file.path(.wd,'out/single_species_models/area_trait')
  .traitPF <- file.path(.wd, 'raw_data/anthropause_data_sheet.csv')
  
  .cores <- 4
  .minsp <- 10
  .iter <- 500
  .thin <- 1
  
} else {
  library(docopt)
  
  ag <- docopt(doc, version = '0.1\n')
  
  .wd <- getwd()
  .cores <- ag$cores
  .minsp <- ag$minsp
  .iter  <- ag$iter
  .thin <- ag$thin
  
  source(file.path(.wd,'analysis/src/funs/input_parse.r'))
  
  .datPF <- makePath(ag$dat)
  .outP <- makePath(ag$out)
  .traitPF <- makePath(ag$trait)
  
}

# ==== Setup ====

#---- Initialize Environment ----#

t0 <- Sys.time()

source(file.path(.wd,'analysis/src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    # library(iterators)
    library(brms)
    library(tidyverse)
    library(lubridate)
    library(glue)
    library(foreach)
    library(doMC)
    
  }))

# Manage package conflicts
conflict_prefer("accumulate", "foreach")
conflict_prefer("ar", "brms")
conflict_prefer("lag", "stats")
conflict_prefer("when", "foreach")

# load breezy functions
source(file.path(.wd,'analysis/src/funs/auto/breezy_funs.r'))

# check arg inputs
.minsp <- ifelse(is.null(.minsp), 10, as.numeric(.minsp))
.iter <- ifelse(is.null(.iter), 3000, as.numeric(.iter))
.thin <- ifelse(is.null(.thin), 4, as.numeric(.thin))

#---- Load data ----#
message("Loading data...")

# load and process data
size <- read_csv("out/dbbmm_size.csv") %>%
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>%
  distinct()

message("Loading trait data...")
traits <- read_csv(.traitPF) %>% 
  mutate(mig_code = case_when(migratory == "Partial" ~ "Partial Migrant",
                              migratory == "Complete" ~ "Complete Migrant",
                              migratory == "migratory" ~ "Complete Migrant",
                              migratory == "non-migratory" ~ "Resident"))
# get ind count per species
sp_sum <- size %>%
  group_by(species) %>%
  summarize(nind = length(unique(ind_id))) %>%
  filter(nind > .minsp) %>% #require a minimum of 10 individuals
  pull(species)

# combine data
dat <- size %>% 
  # filter(species %in% sp_sum) %>% 
  left_join(traits, by = c("species" = "Species")) %>% 
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
         study_f = as.factor(study_id)) # make study id factor)


mammals <- dat %>% 
  filter(Family == "Cervidae" | Family == "Canidae" | Family == "Ursidae" |
           Family == "Felidae" | Family == "Antilocapridae")

birds <- dat %>% 
  filter(Family == "Accipitridae" | Family == "Falconidae" | Family == "Gruidae" |
           Family == "Cathartidae" | Family == "Anatidae" | Family == "Ardeidae" | 
           Family == "Corvidae" | Family == "Rallidae") 

# ==== Perform analysis ====



#-- Mobility by Class --#

form_sg <-  bf(log_area ~ sg_scale*class + ndvi_scale + lst_scale + (1 |species/grp) + ar(time = wk, gr = grp))
message("Fitting models with formula:")
print(form)

message("Starting model...")

# fit model
mod_sg <- brm(
  form_sg,
  data = dat,
  # family = Gamma(link = "log"),
  inits = 0,
  cores = .cores,
  iter = .iter,
  thin = .thin
)

message("Saving model...")
#stash results into named list
out_sg <- list(
  data = dat,
  model = mod_sg
)

#write out results
save(out_sg, file = glue("{.outP}/bird-v-mammal_sg_mod_{Sys.Date()}.rdata"))



#-- Modification by Class --#

form_ghm <-  bf(log_area ~ ghm_scale*class + ndvi_scale + lst_scale + (1 |species/grp) + ar(time = wk, gr = grp))
message("Fitting models with formula:")
print(form)

message("Starting model...")

# fit model
mod_ghm <- brm(
  form_ghm,
  data = dat,
  # family = Gamma(link = "log"),
  inits = 0,
  cores = .cores,
  iter = .iter,
  thin = .thin
)

message("Saving model...")
#stash results into named list
out_ghm <- list(
  data = dat,
  model = mod_ghm
)

#write out results
save(out_ghm, file = glue("{.outP}/bird-v-mammal_ghm_mod_{Sys.Date()}.rdata"))


#-- SG + GHM by Class --#

form_sg_ghm <-  bf(log_area ~ (sg_scale + ghm_scale)*class + ndvi_scale + lst_scale + (1 |species/grp) + ar(time = wk, gr = grp))
message("Fitting models with formula:")
print(form)

message("Starting model...")

# fit model
mod_sg_ghm <- brm(
  form,
  data = dat,
  # family = Gamma(link = "log"),
  inits = 0,
  cores = .cores,
  iter = .iter,
  thin = .thin
)

message("Saving model...")
#stash results into named list
out_sg_ghm <- list(
  data = dat,
  model = mod_sg_ghm
)

#write out results
save(out, file = glue("{.outP}/bird-v-mammal_sg-ghm_mod_{Sys.Date()}.rdata"))

# ==== Finalize script ====
message(glue('Script complete in {diffmin(t0)} minutes'))


