suppressWarnings(
  suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(ggthemes)
    library(patchwork)
    library(brms)
    library(grid)
    library(cowplot)
  }))

# col <- "#7552A3"
col <- "#D3D3D3"
# pal <- c("#7552A3", "#E66100")
pal <- c("#7552A3", "#CEBEDA")


#---- Load and name data ----#

load("out/single_species_models/area_additive/Puma concolor_2022-09-02.rdata")
area_add <- out

load("out/single_species_models/niche_additive/Puma concolor_2022-09-02.rdata")
niche_add <- out

load("out/single_species_models/area/Puma concolor_2022-09-01.rdata")
area_int <- out

load("out/single_species_models/niche/Puma concolor_2022-09-02.rdata")
niche_int <- out




#---- Area Plots ----#


#-- Mobility --#

# Get conditional effects
area_sg <- conditional_effects(x=area_add$model,
                               effects = "sg_norm",
                               re_formula = NA)

(area_sg_ce_plot <-  plot(area_sg, plot = F,
                          line_args = list("se" = T,
                                           "color" = col,
                                           "fill" = col))[[1]] +
    # scale_color_manual(values = palnew[3])+
    # theme_minimal() +
    xlab("Mobility") +
    ylab("Space Use")+
    # ggtitle("Puma concolor")+
    theme_cowplot()  +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          aspect.ratio = 1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    # labs(x  = "Day of year", y = "Mobility") +
    NULL
)
ggsave(plot = area_sg_ce_plot, filename = "figs/cougars/cougar_area_add_sg.png",
       width = 4, height = 4, bg = "white")

#-- Modification --#

# Get conditional effects
area_ghm <- conditional_effects(x=area_add$model,
                                effects = "ghm_scale",
                                re_formula = NA)

(area_ghm_ce_plot <-  plot(area_ghm, plot = F,
                           line_args = list("se" = T,
                                            "color" = col,
                                            "fill" = col))[[1]] +
    # scale_color_manual(values = palnew[3])+
    # theme_minimal() +
    xlab("Modification") +
    ylab("Space Use")+
    # ggtitle("Puma concolor")+
    theme_cowplot()  +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          aspect.ratio = 1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    # labs(x  = "Day of year", y = "Mobility") +
    NULL
)
ggsave(plot = area_ghm_ce_plot, filename = "figs/cougars/cougar_area_add_ghm.png",
       width = 4, height = 4, bg = "white")


#-- Interaction --#

## get observed quantiles of ghm to set "low" and "high" human mod
area_ghmq <- quantile(area_int$data$ghm_scale, probs = c(0.10, 0.90), na.rm = T)

# Conditional Effects Plot for interaction
area_ce_int <- conditional_effects(x=area_int$model, 
                                   effects = "sg_norm:ghm_scale",
                                   int_conditions = list(ghm_scale = area_ghmq),
                                   re_formula = NA)

(area_int_ce_plot <-  plot(area_ce_int, plot = FALSE,
                           line_args = list("se"=T,
                                            "alpha" = 0.2))[[1]] +
    scale_color_manual(values = pal, name = "Modification",
                       labels = c("High", "Low")) +
    scale_fill_manual(values = pal, name = "Modification",
                      labels = c("High", "Low")) +
    xlab("Mobility") +
    ylab("Niche Breadth")+
    theme_cowplot()  +
    theme(legend.position = "none",
          legend.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          aspect.ratio = 1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    NULL
)
ggsave(plot = area_int_ce_plot, filename = "figs/cougars/cougar_area_int.png",
       width = 4, height = 4, bg = "white")



#---- Niche Plots ----#


#-- Mobility --#

# Get conditional effects
niche_sg <- conditional_effects(x=niche_add$model,
                                effects = "sg_norm",
                                re_formula = NA)

(niche_sg_ce_plot <-  plot(niche_sg, plot = F,
                           line_args = list("se" = T,
                                            "color" = col,
                                            "fill" = col))[[1]] +
    # scale_color_manual(values = palnew[3])+
    # theme_minimal() +
    xlab("Mobility") +
    ylab("Niche Breadth")+
    # ggtitle("Puma concolor")+
    theme_cowplot()  +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          aspect.ratio = 1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    # labs(x  = "Day of year", y = "Mobility") +
    NULL
)
ggsave(plot = niche_sg_ce_plot, filename = "figs/cougars/cougar_niche_add_sg.png",
       width = 4, height = 4, bg = "white")

#-- Modification --#

# Get conditional effects
niche_ghm <- conditional_effects(x=niche_add$model,
                                 effects = "ghm_scale",
                                 re_formula = NA)

(niche_ghm_ce_plot <-  plot(niche_ghm, plot = F,
                            line_args = list("se" = T,
                                             "color" = col,
                                             "fill" = col))[[1]] +
    # scale_color_manual(values = palnew[3])+
    # theme_minimal() +
    xlab("Modification") +
    ylab("Niche Breadth")+
    # ggtitle("Puma concolor")+
    theme_cowplot()  +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          aspect.ratio = 1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    # labs(x  = "Day of year", y = "Mobility") +
    NULL
)
ggsave(plot = niche_ghm_ce_plot, filename = "figs/cougars/cougar_niche_add_ghm.png",
       width = 4, height = 4, bg = "white")


#-- Interaction --#

## get observed quantiles of ghm to set "low" and "high" human mod
niche_ghmq <- quantile(niche_int$data$ghm_scale, probs = c(0.10, 0.90), na.rm = T)

# Conditional Effects Plot for interaction
niche_ce_int <- conditional_effects(x=niche_int$model, 
                                    effects = "sg_norm:ghm_scale",
                                    int_conditions = list(ghm_scale = niche_ghmq),
                                    re_formula = NA)

(niche_int_ce_plot <-  plot(niche_ce_int, plot = FALSE,
                            line_args = list("se"=T,
                                             "alpha" = 0.2))[[1]] +
    scale_color_manual(values = pal, name = "Modification",
                       labels = c("High", "Low")) +
    scale_fill_manual(values = pal, name = "Modification",
                      labels = c("High", "Low")) +
    xlab("Mobility") +
    ylab("Niche Breadth")+
    theme_cowplot()  +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          aspect.ratio = 1) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    NULL
)
ggsave(plot = niche_int_ce_plot, filename = "figs/cougars/cougar_niche_int.png",
       width = 4, height = 4, bg = "white")
