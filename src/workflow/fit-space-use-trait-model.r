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
  .outP <- file.path(.wd,'out/single_species_models/area')
  .traitPF <- file.path(.wd, 'raw_data/anthropause_data_sheet.csv')
  
  .cores <- 4
  .minsp <- 10
  .iter <- 5000
  .thin <- 4
  
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
  mutate(
    log_area = log(area), #get log of weekly area use
    log_area_scale = scale(log_area), # standardize it
    sg_norm = scale(sg / cbg_area), # normalize safegraph data by size of the CBG
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
    study_f = as.factor(study_id) # make study id factor
  ) %>%
  distinct()

message("Loading trait data...")
traits <- read_csv(.traitPF)

# combine data
dat <- size %>% 
  left_join(traits, by = c("species" = "Species")) %>% 
  mutate(diet = case_when(Diet.PlantO >= 50 ~ "Herbivore",
                          Diet.Fruit >= 50 ~ "Frugivore",
                          Diet.Scav >= 50 ~ "Savenger",
                          Diet.Inv >= 50 ~ "Insectivore",
                          (Diet.Vend+Diet.Vect+Diet.Vfish) >= 50 ~ "Carnivore",
                          TRUE ~ "Omnivore"))


# ==== Perform analysis ====

#declare model form
form <-  bf(log_area_scale ~ sg_norm*ghm_scale + ndvi_scale + lst_scale + (1 |species/grp) + ar(time = wk, gr = grp))
message("Fitting models with formula:")
print(form)

message("Strating model...")

# fit model
mod <- brm(
  form,
  data = dat,
  # family = Gamma(link = "log"),
  inits = 0,
  cores = .cores,
  iter = .iter,
  thin = .thin
)

#stash results into named list
out <- list(
  data = dat,
  model = mod
)

#write out results
save(out, file = glue("{.outP}/size_trait_mod_{Sys.Date()}.rdata"))


# ==== Finalize script ====
message(glue('Script complete in {diffmin(t0)} minutes'))


