
# Same set up as the fit-space-use-interactive-models.r script. 
# (to assess which data make it into there and what makes it through the loop)
if(interactive()) {
  
  .wd <- getwd()
  
  .datPF <- file.path(.wd,'out/dbbmm_size.csv')
  .outP <- file.path(.wd,'out/single_species_models/area_interactive')
  
  .cores <- 20
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
  
}

# ==== Setup ====

#---- Initialize Environment ----#

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
conflict_prefer("expand", "tidyr")
conflict_prefer("has_name", "tibble")
conflict_prefer("pack", "tidyr")
conflict_prefer("unpack", "tidyr")

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
  mutate(ind_f = as.factor(ind_id))%>%  # create factor version of ind for REs)
  mutate(species = case_when( # correct species names
    study_id == 1442516400 ~ "Anser caerulescens",
    study_id == 1233029719 ~ "Odocoileus virginianus",
    study_id == 1631574074 ~ "Ursus americanus",
    study_id == 1418296656 ~ "Numenius americanus",
    study_id == 474651680  ~ "Odocoileus virginianus",
    study_id == 1044238185 ~ "Alces alces",
    study_id == 2548691779 ~ "Odocoileus hemionus", # !!! CORRECT IN WORKFLOW !!! #
    study_id == 2575515057 ~ "Cervus elaphus",  # !!! CORRECT IN WORKFLOW !!! #
    TRUE ~ species
  ))%>% 
  mutate(species = case_when(
    species == "Chen caerulescens" ~ "Anser caerulescens",
    TRUE ~ species
  )) %>% 
  distinct()

# #investiagte NA species
# size %>% 
#   filter(is.na(species)) %>% 
#   pull(study_id) %>% 
#   unique()

# get ind count per species
sp_sum <- size %>%
  group_by(species) %>%
  summarize(nind = length(unique(ind_f))) %>%
  filter(nind > .minsp) #require a minimum of 10 individuals
sp_sum$species

# Filter complete cases, then species with too few inds
print(complete <- size %>% 
        select(species, ind_id, tmax, ndvi, area, sg, ghm) %>% 
        filter(complete.cases(.)) %>% 
        group_by(species) %>%
        summarize(nind = length(unique(ind_id)))%>%
        filter(nind > .minsp), n= 100)
complete$species

(failed_size_spp_ls <- sp_sum %>% 
    filter(species %notin% complete$species) %>% 
    pull(species))

# filter to only compleet data going in
complete_dat <- size %>% 
  select(species, ind_id, tmax, ndvi, area, sg, ghm) %>% 
  filter(complete.cases(.))

failed_size_dat <- size %>% 
  filter(species %in% failed_size_spp_ls)


####----    Check on Niche Data.   ----#####


# load and process data
breadth <- read_csv('out/niche_determinant_anthropause.csv') %>%
  mutate(scientificname = case_when( # correct species names
    studyid == 1442516400 ~ "Anser caerulescens",
    studyid == 1233029719 ~ "Odocoileus virginianus",
    studyid == 1631574074 ~ "Ursus americanus",
    studyid == 1418296656 ~ "Numenius americanus",
    studyid == 474651680  ~ "Odocoileus virginianus",
    studyid == 1044238185 ~ "Alces alces",    
    studyid == 2548691779 ~ "Odocoileus hemionus", # !!! CORRECT IN WORKFLOW !!! #
    studyid == 2575515057 ~ "Cervus elaphus",  # !!! CORRECT IN WORKFLOW !!! #
    TRUE ~ scientificname
  ))%>% 
  mutate(scientificname = case_when(
    scientificname == "Chen caerulescens" ~ "Anser caerulescens",
    TRUE ~ scientificname
  )) %>% 
  distinct() %>%
  mutate(tmax_mvnh = tmax,
         # tmin_mvnh = tmin,
         elev_mvnh = elev,
         ndvi_mvnh = ndvi) %>%
  select(!c(tmax, elev, ndvi)) %>%
  filter(studyid != 351564596) %>%
  filter(studyid != 1891587670) %>%
  mutate(ind_f = as.factor(individual)) %>%  # create factor version of ind for REs)
  
  left_join(size, by = c("scientificname" = "species", 
                           "individual" = "ind_id", 
                           "year" = "year", 
                           "studyid" = "study_id", 
                           "week" = "wk",
                         "ind_f" = "ind_f"))

# find and fix the NAs
breadth %>% 
  filter(is.na(scientificname)) %>% 
  pull(studyid) %>% 
  unique()

# get ind count per species
nich_sp_sum <- breadth %>%
  group_by(scientificname) %>%
  summarize(nind = length(unique(ind_f))) %>%
  filter(nind > .minsp) #require a minimum of 10 individuals
print(nich_sp_sum, n = 50)

# get list of complete species
niche_complete <- breadth %>% 
  select(scientificname, ind_f, tmax, ndvi, area, sg, ghm) %>% 
  filter(complete.cases(.)) %>% 
  group_by(scientificname) %>%
  summarize(nind = length(unique(ind_f)))%>%
  filter(nind > .minsp)
print(niche_complete, n = 50)


comb_complete <- rbind(niche_complete %>% rename(species = scientificname), complete)

unique(comb_complete$species)
