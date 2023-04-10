#!/usr/bin/env Rscript --vanilla

# ==== Input parameters ====

'
Usage:
fit-space-use-models.r <nichedat> <vardat> <trait> <out> <cores> [<minsp> <iter> <thin>]
fit-space-use-models.r (-h | --help)


Parameters:
nichedat: path to input niche breadth csv file.
vardat: path to dbbmm data csv with covars
out: path to output directory.
cores: number of HPC cores 
iter: MCMC iterations
thin: MCMC thin rate

Options:
-h --help           Show this screen.
-v --version        Show version.

' -> doc

if(interactive()) {
  
  .wd <- getwd()
  
  .datPF <- file.path(.wd,'out/niche_determinant_anthropause.csv')
  .varPF <- file.path(.wd, 'out/dbbmm_size.csv')
  .traitPF <- file.path(.wd, 'raw_data/anthropause_data_sheet.csv')
  .outP <- file.path(.wd,'out/single_species_models/area_trait')
  
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

dbbmms <- read_csv(.varPF) %>% 
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
  ))

# load and process data
breadth <- read_csv(.datPF) %>%
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
         tmin_mvnh = tmin,
         lst_mvnh = lst,
         ndvi_mvnh = ndvi) %>%
  select(!c(tmax, tmin, lst, ndvi)) %>%
  filter(studyid != 351564596) %>%
  filter(studyid != 1891587670) %>%
  mutate(ind_f = as.factor(individual)) %>%  # create factor version of ind for REs)
  left_join(dbbmms, by = c("scientificname" = "species", 
                           "individual" = "ind_id", 
                           "year" = "year", 
                           "studyid" = "study_id", 
                           "week" = "wk")) %>% 
    mutate(sqrt_breadth = sqrt(total), #get log of weekly area use
           sqrt_breadth_scale = scale(sqrt_breadth), # standardize it
           sg_norm = scale(sg / cbg_area), # normalize safegraph data by size of the CBG
           # log_sg_norm = log(sg_norm),
           ghm_scale = scale(ghm),
           ndvi_scale = scale(ndvi.y),
           lst_scale = scale(lst.y),
           
           grp = paste(ind_f, year, sep = "_"), # create indXyr grouping factor
           # trt_new = gsub('_.*','',trt),
           year_f = factor(year), # create year factor
           # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
           # sp2 = gsub(" ", "_", species),
           wk_n = as.numeric(substring(week, 2)), # extract week number
           ts = parse_date_time(paste(year, week, 01, sep = "-"), "%Y-%U-%u"), # make better date format
           study_f = as.factor(studyid)) # make study id factor)

# get ind count per species
sp_sum <- breadth %>%
  group_by(scientificname) %>%
  summarize(nind = length(unique(ind_f))) %>%
  filter(nind > .minsp) #require a minimum of 10 individuals

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
mammals <- dat %>% 
  filter(Family == "Cervidae" | Family == "Canidae" | Family == "Ursidae" |
           Family == "Felidae" | Family == "Antilocapridae")

birds <- dat %>% 
  filter(Family == "Accipitridae" | Family == "Falconidae" | Family == "Gruidae" |
           Family == "Cathartidae" | Family == "Anatidae" | Family == "Ardeidae" | 
           Family == "Corvidae" | Family == "Rallidae") 

# ==== Perform analysis ====

# declare model form
form <-  bf(log_area_scale ~ sg_norm*ghm_scale*diet + ndvi_scale + lst_scale + (1 |species/grp) + ar(time = wk, gr = grp))
message("Fitting models with formula:")
print(form)

message("Strating model...")

# fit model
mod_mammal <- brm(
  form,
  data = mammals,
  # family = Gamma(link = "log"),
  inits = 0,
  cores = .cores/2,
  iter = .iter,
  thin = .thin
)

#stash results into named list
out_mammals <- list(
  data = dat,
  model = mod
)

#write out results
save(out, file = glue("{.outP}/mammal_size_trait_mod_{Sys.Date()}.rdata"))

# fit model
mod_birds <- brm(
  form,
  data = birds,
  # family = Gamma(link = "log"),
  inits = 0,
  cores = .cores/2,
  iter = .iter,
  thin = .thin
)

#stash results into named list
out_birds <- list(
  data = dat,
  model = mod
)

#write out results
save(out, file = glue("{.outP}/bird_size_trait_mod_{Sys.Date()}.rdata"))

# ==== Finalize script ====
message(glue('Script complete in {diffmin(t0)} minutes'))


