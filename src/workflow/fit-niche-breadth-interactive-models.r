#!/usr/bin/env Rscript 
#
# ==== Input parameters ====

'
Usage:
fit-niche-breadth-models.r <nichedat> <vardat> <out> <cores> [<minsp> <iter> <thin>]
fit-niche-breadth-models.r (-h | --help)


Parameters:
nichedat: path to input niche breadth csv file.
vardat: path to dbbmm data csv with covars
out: path to output directory.
cores: the number of HPC cores to request

Options:
-h --help           Show this screen.
-v --version        Show version.

' -> doc

if(interactive()) {
  
  .wd <- getwd()
  
  .datPF <- file.path(.wd,'out/niche_determinant_anthropause.csv')
  .varPF <- file.path(.wd, "out/dbbmm_size.csv")
  .outP <- file.path(.wd,'out/single_species_models/niche')
  
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
  
  source(file.path(.wd,'src/funs/input_parse.r'))
  
  .datPF <- makePath(ag$nichedat)
  .varPF <- makePath(ag$vardat)
  .outP <- makePath(ag$out)
  
}

# ==== Setup ====

#---- Initialize Environment ----#

t0 <- Sys.time()

source(file.path(.wd,'src/startup.r'))

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
options(scipen = 999)
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
source(file.path(.wd,'src/funs/auto/breezy_funs.r'))

# check arg inputs
.minsp <- ifelse(is.null(.minsp), 10, as.numeric(.minsp))
.iter <- ifelse(is.null(.iter), 3000, as.numeric(.iter))
.thin <- ifelse(is.null(.thin), 4, as.numeric(.thin))

#---- Load data ----#
message("Loading data...")

dbbmms <- read_csv(.varPF) %>% 
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
  ))

# load and process data
breadth <- read_csv(.datPF) %>%
  mutate(scientificname = case_when( # correct species names
    studyid == 1442516400 ~ "Anser caerulescens",
    studyid == 1631574074 ~ "Ursus americanus",
    studyid == 1418296656 ~ "Numenius americanus",
    studyid == 474651680  ~ "Odocoileus virginianus",
    studyid == 1044238185 ~ "Alces alces",
    TRUE ~ scientificname
  ))%>% 
  mutate(scientificname = case_when(
    scientificname == "Chen caerulescens" ~ "Anser caerulescens",
    scientificname == "Chen rossii" ~ "Anser rossii",
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
  left_join(dbbmms, by = c("scientificname" = "species", 
                           "individual" = "ind_id", 
                           "year" = "year", 
                           "studyid" = "study_id", 
                           "week" = "wk")) %>%
  # exclude some individuals per data owner's request
  filter(individual != 3418130234, 
        individual != 3418130244, 
        individual != 3418130245, 
        individual != 3418130249)

# get ind count per species
sp_sum <- breadth %>%
  group_by(scientificname) %>%
  summarize(nind = length(unique(ind_f))) %>%
  filter(nind >= .minsp) #require a minimum of 3 individuals

# ==== Start cluster and register backend ====
registerDoMC(.cores)

# ==== Perform analysis ====

#declare model form
form <-  bf(breadth_scale ~ log_area_scale + sg_norm*ghm_scale + ndvi_scale + tmax_scale + (1 |grp) + ar(time = week, gr = grp))
# form <-  bf(breadth_scale ~ sg_norm*ghm_scale + (1 |grp) + ar(time = week, gr = grp))
message("Fitting models with formula:")
print(form)

# loop through species
foreach(i = 1:nrow(sp_sum), .errorhandling = "pass", .inorder = F) %dopar% {
  
  # get focal species
  sp <- sp_sum$scientificname[i]
  
  
  message(glue("Strating model for {sp}..."))
  
  dat <- breadth %>%
    filter(scientificname == sp,
           !is.infinite(total)) %>% 
    mutate(
      log_area = log(area), #get log of weekly area use
      log_area_scale = scale(log_area), # standardize it
      # sqrt_breadth = sqrt(total), #get sqrt weekly niche breadth
      breadth_scale = scale(log(total+0.000000000000000000001)), # standardize it
      sg_norm = scale(sg / cbg_area), # normalize safegraph data by size of the CBG
      # log_sg_norm = log(sg_norm),
      ghm_scale = scale(ghm),
      ndvi_scale = scale(ndvi),
      # lst_scale = scale(lst),
      tmax_scale = scale(tmax),
      # tmin_scale = scale(tmin),
      grp = paste(ind_f, year, sep = "_"), # create indXyr grouping factor
      # trt_new = gsub('_.*','',trt),
      year_f = factor(year), # create year factor
      # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
      # sp2 = gsub(" ", "_", species),
      wk_n = as.numeric(substring(week, 2)), # extract week number
      ts = parse_date_time(paste(year, week, 01, sep = "-"), "%Y-%U-%u"), # make better date format
      study_f = as.factor(studyid) # make study id factor
    ) %>%
    distinct()    
  
  # fit model
  mod <- brm(
    form,
    data = dat,
    # family = Gamma(link = "log"),
    init = 0,
    cores = 1,
    iter = .iter,
    thin = .thin,
    control = list(adapt_delta = 0.95)
  )
  
  #stash results into named list
  out <- list(
    species = sp,
    data = dat,
    model = mod
  )
  
  #write out results
  save(out, file = glue("{.outP}/{sp}_{Sys.Date()}.rdata"))
  
} # i

# ==== Finalize script ====
message(glue('Script complete in {diffmin(t0)} minutes'))


