library(meta)
library(tidyverse)
library(brms)
library(parameters)
library(bayestestR)

# conflicts_prefer(brms::ar)
# conflicts_prefer(tibble::has_name)
# conflicts_prefer(dplyr::lag)

# declare vector of birds
birds <- c("Anser caerulescens", "Aquila chrysaetos", "Ardea alba", "Corvus corax", 
           "Grus canadensis","Haliaeetus leucocephalus", "Numenius americanus")

sg_dat <- read_csv("out/niche_sg_marginal_2023-09-27.csv") %>% 
  mutate(se = uncertainty/1.96,
         wei = 1/uncertainty^2,
         wei_norm = wei/sum(wei),
         class = case_when(species %in% birds ~ "bird",
                           TRUE ~ "mammal")) 

brm_sg <- brm(
  Estimate | weights(wei, scale = TRUE) ~ 0 + class, 
  data = sg_dat, 
  cores = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99)
)

brm_sg
plot(brm_sg)
pd_sg <- p_direction(brm_sg) %>% 
  mutate(Parameter = case_when(Parameter == "b_classbird" ~ "classbird",
                               Parameter == "b_classmammal" ~ "classmammal"))
out_sg <- fixef(brm_sg) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "Parameter") %>% 
  left_join(pd_sg)
write_csv(out_sg, "out/niche_meta_sg.csv")

# GHM
ghm_dat <- read_csv("out/niche_ghm_marginal_2023-09-27.csv") %>% 
  mutate(se = uncertainty/1.96,
         wei = 1/uncertainty^2,
         wei_norm = wei/sum(wei),
         class = case_when(species %in% birds ~ "bird",
                           TRUE ~ "mammal")) 



brm_ghm <- brm(
  Estimate | weights(wei, scale = TRUE) ~ 0 + class, 
  data = ghm_dat, 
  cores = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99)
)

brm_ghm
plot(brm_ghm)
pd_ghm <- p_direction(brm_ghm) %>% 
  mutate(Parameter = case_when(Parameter == "b_classbird" ~ "classbird",
                               Parameter == "b_classmammal" ~ "classmammal"))
out_ghm <- fixef(brm_ghm) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "Parameter") %>% 
  left_join(pd_ghm)
write_csv(out_ghm, "out/niche_meta_ghm.csv")
