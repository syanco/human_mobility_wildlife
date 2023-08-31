library(meta)
library(tidyverse)
library(brms)
library(parameters)

# conflicts_prefer(brms::ar)
# conflicts_prefer(tibble::has_name)
# conflicts_prefer(dplyr::lag)

# declare vector of birds
birds <- c("Anser caerulescens", "Aquila chrysaetos", "Ardea alba", "Corvus corax", "Haliaeetus leucocephalus")

sg_dat <- read_csv("out/area_sg_marginal_2023-06-21.csv") %>% 
  mutate(se = uncertainty/1.96,
         wei = 1/uncertainty^2,
         wei_norm = wei/sum(wei),
         class = case_when(species %in% birds ~ "bird",
                           TRUE ~ "mammal")) 

brm_sg <- brm(
  Estimate | weights(wei, scale = TRUE) ~ 0 + class + (1 | species), 
  data = sg_dat, 
  cores = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99)
)

brm_sg
plot(brm_sg)
out_sg <- parameters(brm_sg)
write_csv(out_sg, "out/area_meta_sg.csv")

# GHM
ghm_dat <- read_csv("out/area_ghm_marginal_2023-06-21.csv") %>% 
  mutate(se = uncertainty/1.96,
        wei = 1/uncertainty^2,
        wei_norm = wei/sum(wei),
        class = case_when(species %in% birds ~ "bird",
                          TRUE ~ "mammal")) 



brm_ghm <- brm(
  Estimate | weights(wei, scale = TRUE) ~ 0 + class + (1 | species), 
  data = ghm_dat, 
  cores = 4,
  iter = 10000,
  control = list(adapt_delta = 0.99)
)

brm_ghm
plot(brm_ghm)
out_ghm <- parameters(brm_ghm)
write_csv(out_ghm, "out/area_meta_ghm.csv")
