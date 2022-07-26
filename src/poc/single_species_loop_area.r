# quick script to loop through species and fit a model for each...

library(brms)
library(tidyverse)
library(lubridate)
library(glue)
library(foreach)
library(doMC)


registerDoMC(10)

size <- read_csv("out/dbbmm_size.csv") %>%
  mutate(
    # trt_new = gsub('_.*','',trt),
    year_f = factor(year),
    # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
    sp2 = gsub(" ", "_", species),  
    ind_f = as.factor(ind_id),
    wk_n = as.numeric(substring(wk, 2)),
    ts = parse_date_time(paste(year, wk, 01, sep = "-"), "%Y-%U-%u")
  ) %>%
  distinct()

#- Load ain the species trait data
traits <- read_csv("raw_data/anthropause_data_sheet.csv") %>% 
  mutate(mig_mod = case_when(migratory == "Partial" ~ "migratory",
                             .$migratory == "non-migratory" ~ "non-migratory",
                             .$migratory == "Complete" ~ "migratory",
                             .$migratory == "migratory" ~ "migratory",
                             .$migratory == "Semi-nomadic" ~ "non-migratory",
                             .$migratory == "Unknown" ~ NA_character_))

size <- size %>% 
  left_join(traits, by = c("species" = "Species"))


size_res <- size %>%
  # filter(mig_mod == "non-migratory") %>%
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>%
  mutate(
    log_area = log(area),
    sg_norm = area / cbg_area,
    log_sg_norm = log(sg_norm),
    grp = paste(ind_f, year, sep = "_")
  )

# get ind count per species
sp_sum <- size_res %>%
  group_by(sp2) %>%
  summarize(nind = length(unique(ind_f))) %>%
  filter(nind >10) #require a minimum of 10 individuals


#declare model form
form <-  bf(log_area ~ sg_norm*ghm + ndvi + lst + (1 |grp) + ar(time = wk, gr = grp))

# loop throuigh species
foreach(i = 1:nrow(sp_sum), .errorhandling = "pass", .inorder = F) %dopar% {
  
  # get focal species
  sp <- sp_sum$sp2[i]
  
  #filter data
  dat <- size_res %>%
    filter(sp2 == sp)
  
  # fit model
  mod <- brm(
    form,
    # data2 = list(phylo_vcov = phylo_vcov),
    data = dat,
    # family = Gamma(link = "log"),
    inits = 0,
    cores = 4,
    iter = 3000,
    thin = 4
  )
  
  #stash results into named list
  out <- list(
    species = sp,
    data = dat,
    model = mod
  )
  
  #write out results
  save(out, file = glue("out/single_species_models/area/{sp}_{Sys.Date()}.rdata"))
}