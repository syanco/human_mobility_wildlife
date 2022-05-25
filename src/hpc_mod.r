library(brms)
library(tidyverse)
library(lubridate)


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
  filter(mig_mod == "non-migratory") %>%
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>%
  mutate(
    log_area = log(area),
    sg_norm = area / cbg_area,
    log_sg_norm = log(sg_norm),
    grp = paste(ind_f, year, sep = "_")
  )

form <-  bf(log_area ~ sg_norm*ghm + ndvi + lst + (1 + sg_norm*ghm + ndvi + lst | sp2) + (1 |sp2:grp) + ar(time = wk , gr = grp))

mod_res_rslopes <- brm(
  form,
  # data2 = list(phylo_vcov = phylo_vcov),
  data = size_res,
  # family = Gamma(link = "log"),
  inits = 0,
  cores = 4,
  iter = 5000,
  thin = 4
)

# pp_check(mod_res)
# plot(mod_res)

save(mod_res_rslopes, file = "out/quick_mod_rslopes_05162022.rdata")
