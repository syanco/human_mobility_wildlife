# quick script to loop through species and fit a model for each...

library(brms)
library(tidyverse)
library(lubridate)
library(glue)
library(foreach)
library(doMC)


registerDoMC(10)

breadth <- read_csv("out/niche_comb.csv") %>%
  mutate(
    # trt_new = gsub('_.*','',trt),
    year_f = factor(year),
    # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
    sp2 = gsub(" ", "_", species),  
    ind_f = as.factor(individual),
    wk_n = as.numeric(substring(week, 2)),
    ts = parse_date_time(paste(year, week, 01, sep = "-"), "%Y-%U-%u")
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

breadth <- breadth %>% 
  left_join(traits, by = c("species" = "Species"))


breadth_res <- breadth %>%
  # filter(mig_mod == "non-migratory") %>%
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>%
  filter(!is.infinite(total),
         !is.na(total)) %>%
  mutate(
    log_area = log(area),
    sg_norm = area / cbg_area,
    log_sg_norm = log(sg_norm),
    grp = paste(ind_f, year, sep = "_"),
    tot_norm = sqrt(total)
  )

# get ind count per species
sp_sum <- breadth_res %>%
  group_by(sp2) %>%
  summarize(nind = length(unique(ind_f))) %>%
  filter(nind >10) #require a minimum of 10 individuals


#decalre model form
form <-  bf(tot_norm ~ sg_norm*ghm + (1 |grp) + ar(time = week, gr = grp))

# loop throuigh species
foreach(i = 1:nrow(sp_sum), .errorhandling = "pass", .inorder = F) %dopar% {
  
  # get focal species
  sp <- sp_sum$sp2[i]
  
  #filter data
  dat <- breadth_res %>%
    filter(sp2 == sp) %>%
    filter(!is.na(total) & !is.na(sg) & !is.na(ghm) & !is.na(ndvi.y), !is.na(lst.y)) 
  
  # fit model
  mod <- brm(
    form,
    # data2 = list(phylo_vcov = phylo_vcov),
    data = dat,
    family = skew_normal(),
    # inits = 0,
    cores = 4,
    iter = 3000,
    thin = 2,
    control = list(adapt_delta = 0.95)
  )
  
  #stash results into named list
  out <- list(
    species = sp,
    data = dat,
    model = mod
  )
  
  #write out results
  save(out, file = glue("out/single_species_models/niche/{sp}_{Sys.Date()}.rdata"))
}
