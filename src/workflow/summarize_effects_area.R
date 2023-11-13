library(meta)
library(tidyverse)
library(brms)
library(parameters)
library(bayestestR)

# conflicts_prefer(brms::ar)
# conflicts_prefer(tibble::has_name)
# conflicts_prefer(dplyr::lag)

# declare vector of birds
birds <- c("Anas acuta", "Anas americana", "Anas clypeata", "Anas crecca",
          "Anas cyanoptera", "Anas platyrhynchos", "Anas strepera", 
          "Anser albifrons", "Anser caerulescens", "Aquila chrysaetos", 
          "Ardea alba", "Chen rossii", "Circus cyaneus", "Corvus corax", 
          "Haliaeetus leucocephalus", "Rallus longirostris")

sg_dat <- read_csv("out/area_sg_marginal_2023-10-17.csv") %>% 
  mutate(se = uncertainty/1.96,
         wei = 1/uncertainty^2,
         wei_norm = wei/sum(wei),
         class = case_when(species %in% birds ~ "bird",
                           TRUE ~ "mammal")) 

sg_dat %>% 
  ggplot(aes(x=Estimate, color = class))+geom_density()

brm_sg <- brm(
  Estimate | weights(wei, scale = TRUE) ~ 0 + class, 
  data = sg_dat ,
  cores = 4,
  iter = 1000000,
  thin = 5,
  init = 0,
  control = list(adapt_delta = 0.99)
)

brm_sg

beepr::beep()
# pp_check(brm_sg)
# plot(brm_sg)
# bayesplot::mcmc_pairs(brm_sg)
(pd_sg <- p_direction(brm_sg) %>% 
  mutate(Parameter = case_when(Parameter == "b_classbird" ~ "classbird",
                   Parameter == "b_classmammal" ~ "classmammal")))
out_sg <- fixef(brm_sg) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "Parameter") %>% 
  left_join(pd_sg)

write_csv(out_sg, "out/area_meta_sg.csv")

# GHM
ghm_dat <- read_csv("out/area_ghm_marginal_2023-10-17.csv") %>% 
  mutate(se = uncertainty/1.96,
        wei = 1/uncertainty^2,
        wei_norm = wei/sum(wei),
        class = case_when(species %in% birds ~ "bird",
                          TRUE ~ "mammal")) 



brm_ghm <- brm(
  Estimate | weights(wei, scale = TRUE) ~ 0 + class + (1|species), 
  data = ghm_dat, 
  cores = 4,
  iter = 100000,
  thin = 3,
  init = 0
  # control = list(adapt_delta = 0.99)
)
beepr::beep()

# plot(brm_ghm)
brm_ghm
(pd_ghm <- p_direction(brm_ghm) %>% 
  mutate(Parameter = case_when(Parameter == "b_classbird" ~ "classbird",
                               Parameter == "b_classmammal" ~ "classmammal")))
out_ghm <- fixef(brm_ghm) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "Parameter") %>% 
  left_join(pd_ghm)
write_csv(out_ghm, "out/area_meta_ghm.csv")
