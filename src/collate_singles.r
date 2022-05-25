# quick script to pull info from single speceis use area models 

library(ggplot2)
library(tidyverse)
library(brms)
library(tidybayes)    # Manipulate Stan objects in a tidy way
library(broom)        # Convert model objects to data frames
library(broom.mixed)  # Convert brms model objects to data frames
library(emmeans)
library(lubridate)
# library(doMC)
# library(foreach)
library(patchwork)
library(glue)
library(ggtheme)
library(magrittr)
library(gridExtra)
library(parameters)

## INITS ##
# .nc <- 4
# registerDoMC(.nc)

fl <- list.files("out/single_species_models/area/")

for(i in 1:length(fl)){
  
  #load the model, species, and data
  load(glue("out/single_species_models/area/{fl[i]}")) 
  
  # produce model summary df
  sum_tab <- out$model %>% 
    parameters(effects = "all", component = "all") %>% as_tibble() %>%
    tableGrob(theme = ttheme_minimal(base_size = 4))
  
  sum_tab_p <- ggplot() + 
    labs(title = "Model Summary") + 
    annotation_custom(sum_tab) + 
    
    theme_minimal()
  
  # produce data summary df
  dat_sum_tab <- out$data %>%
    summarize(
      "No. Individuals" = length(unique(ind_f)),
      "No. Individual-Years" = length(unique(grp)),
      "No. Weeks" = n()) %>%
    pivot_longer(cols = everything(),
                 names_to = "Summary",
                 values_to = "Count") %>%
    tableGrob(theme = ttheme_minimal(base_size = 12))
  
  dat_sum_tab_p <- ggplot() + 
    annotation_custom(dat_sum_tab) + 
    labs(title = "Data Summary") + 
    theme_minimal()
  
  # posterior predictive plot
  ppplot <- pp_check(out$model) +
    ggtitle("Posterior Predictive Check")
  
  # coefficient interval plot
  int_plot <- mcmc_plot(out$model, type = 'intervals') +
    ggtitle("Coefficient Estimates")
  
  # Conditional Effects Plot for interaction
  ce_int <- conditional_effects(x=out$model, 
                                effects = "sg_norm:ghm",
                                # int_conditions = list(ghm = seq(0, 1, by = 0.25)),
                                # rug = T,
                                # points = T,
                                re_formula = NA)
  ce_int_plot <-  plot(ce_int, plot = FALSE)[[1]] +
    theme_minimal() +
    labs(title = "Conditional Effects:",
         subtitle = "Human Mobility x Human Modification") +
    geom_point(data = out$data, aes(x = sg_norm, y = log_area, size = ghm), alpha = 0.2, inherit.aes = F)  
  
  # Conditional Effects Plot for SG
  ce_sg <- conditional_effects(x=out$model, 
                               effects = "sg_norm",
                               # rug = T,
                               # points = T,
                               re_formula = NA)
  ce_sg_plot <-  plot(ce_sg, plot = FALSE)[[1]] +
    theme_minimal() +
    labs(title = "Conditional Effects:",
         subtitle = "Human Mobility") +
    geom_point(data = out$data, aes(x = sg_norm, y = log_area, size = ghm), alpha = 0.2, inherit.aes = F)
  
  
  # Conditional Effects Plot for GHM
  ce_ghm <- conditional_effects(x=out$model, 
                                effects = "ghm",
                                # rug = T,
                                # points = T,
                                re_formula = NA)
  ce_ghm_plot <-  plot(ce_ghm, plot = FALSE)[[1]] +
    theme_minimal() +
    labs(title = "Conditional Effects:",
         subtitle = "Human Modification") +
    geom_point(data = out$data, aes(x = sg_norm, y = log_area, size = sg_norm), alpha = 0.2, inherit.aes = F)
  
  
  
  # get area over time plot
  raw_plot <- ggplot(data = out$data, aes(x = wk, y = log_area, goup = grp, color = year_f)) +
    # geom_point(aes(color = ind_f)) +
    geom_point() +
    geom_line()+
    ggtitle("Space Use Over Time")+
    theme_minimal()
  
  # Put plots together
  out_plot <- (wrap_elements(dat_sum_tab_p) | raw_plot) / 
    ( int_plot | ppplot) / 
    (ce_int_plot) / 
    (ce_sg_plot | ce_ghm_plot) /
    sum_tab_p +
    # (wrap_elements(textGrob("Model Summary"))/wrap_elements(sum_tab)) +
    plot_annotation(title = glue("{out$species}"))& 
    theme(plot.title = element_text(size = 20))
  
  ggsave(filename = glue("out/single_species_models/area_plots/{out$species}.pdf"), 
         plot = out_plot,
         width = 8.5, height = 12.5, units = c("in"))
} #i


#write model summaries ot txt file
sink("out/single_species_models/area_mods.txt", type = "output")
for(i in 1:length(fl)){
  
    
  print(glue("Loading model from {fl[i]}"))
  #load the model, species, and data
  load(glue("out/single_species_models/area/{fl[i]}")) 
  
  print(summary(out$model))
  
  print(glue("\n 
             ############################################################
             \n"))
}
sink()

closeAllConnections()
