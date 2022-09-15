library(brms)
library(tidyverse)
library(tidybayes)

load("out/single_species_models/area_variance_components_mod_2022-09-08.rdata")

traits <- read_csv("raw_data/anthropause_data_sheet.csv")

mod <- out_sg$model

mod
re <- ranef(mod)
fe <- fixef(mod)
re$species
mod$model
re_sp_df <- data.frame(fe) %>% 
  rownames_to_column(var = "param") %>% 
  mutate(est_exp = exp(Estimate),
         species = word(param, 2, sep = "_"),
         species = str_sub(species, 8),
         species = fct_reorder(species, Estimate))

ggplot(re_sp_df) +
  geom_pointinterval(aes(x=Estimate, xmin = Q2.5, xmax = Q97.5, y = species))



re_ind_df <- data.frame(re$`species:grp`) %>% 
  rownames_to_column(var = "grp") %>% 
  mutate(grp = fct_reorder(grp, Estimate.Intercept),
         species = word(grp, 1, sep = "_"),
         ind = fct_reorder(grp, species)) %>% 
  left_join(traits, by = c("species" = "Species"))

ggplot(re_ind_df) +
  geom_pointinterval(aes(x=Estimate.Intercept, xmin = Q2.5.Intercept, xmax = Q97.5.Intercept, y = ind, color = class)) +
  facet_wrap(~species, scales = "free_y") +
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.y=element_blank()) +
  xlab("")+
  NULL



# pp_check(mod)

conditional_effects(mod)

parnames(mod)
mod %>%
  gather_draws() %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = b_cw_Intercept + r_site__cw) %>%
  ungroup() %>%
  mutate(site = str_replace_all(site, "[.]", " ")) %>% 
  
  # plot
  ggplot(aes(x = mu, y = reorder(site, mu))) +
  geom_vline(xintercept = fixef(k_fit_brms)[1, 1], color = "#839496", size = 1) +
  geom_vline(xintercept = fixef(k_fit_brms)[1, 3:4], color = "#839496", linetype = 2) +
  geom_halfeyeh(.width = .5, size = 2/3, fill = "#859900") +
  labs(x = expression("Cottonwood litterfall (g/m^2)"),
       y = "BEMP sites ordered by mean predicted litterfall") +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        