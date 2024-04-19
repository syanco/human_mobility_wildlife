
#-- Init --#

#- Libraries
library(tidyverse)
library(brms)
library(ggthemes)
library(glue)
library(patchwork)

#- Color palette
pal <- c("#F98177","#8895BF")

#-- Load Data --#

#- Niche model
load("out/intra_ind_models/niche_intra_ind_int_mod_2023-09-27.rdata")
niche_out <- out
niche_int_mod <- niche_out$mod

#- Area Model
load("out/intra_ind_models/size_intra_ind_int_mod_2023-09-27.rdata")
area_out <- out
area_int_mod <- out$mod

#-- Plots --#

#- Niche model
niche_ghmq <- quantile(niche_out$data$ghm_diff, probs = c(0.05, 0.95), na.rm = T)

niche_ce_int <- conditional_effects(x=niche_int_mod,
                                    effects = "sg_diff:ghm_diff",
                                    int_conditions = list(ghm_diff = niche_ghmq),
                                    re_formula = NA)
(niche_int_ce_plot <-  plot(niche_ce_int, 
                            plot = F,
                            rug = F,
                            line_args = list("se" = T))[[1]] + 
    scale_color_manual(name ="Changes in human modification",
                       values = rev(pal),
                       labels = c("high",
                                  "low")) +
    scale_fill_manual(name ="Changes in human modification",
                      values = rev(pal),
                      labels = c("high",
                                 "low")) +
    theme_minimal() +
    # xlab(glue("{expression(delta)} Human Mobility")) +
    # xlab(bquote(~Delta~"Human Mobility")) +
    # ylab(bquote(~Delta~"Space Use"))+
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          # legend.position = "none",
          plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 10, 
                                      face = "bold"),
          axis.title.x = element_text(size = 10, 
                                      face = "bold"),
          axis.ticks.x = element_line(color = "#4a4e4d"),
          text = element_text(family = "Arial", color = "#4a4e4d")) +
    labs(x = "Relative change in human mobility", y = "Standardized \n change in niche size", tag = bquote(bold("b")))
)

#- Area Model
area_ghmq <- quantile(area_out$data$ghm_diff, probs = c(0.05, 0.95), na.rm = T)
# area_ghmq <- c(-0.5, 0, 0.5)

area_ce_int <- conditional_effects(x=area_int_mod,
                                   effects = "sg_diff:ghm_diff",
                                   int_conditions = list(ghm_diff = area_ghmq),
                                   re_formula = NA)
(area_int_ce_plot <-  plot(area_ce_int, 
                           plot = F,
                           rug = F,
                           line_args = list("se" = T))[[1]] + 
    scale_color_manual(name ="Changes in human modification",
                       values = rev(pal),
                       labels = c("high",
                                  "low")) +
    scale_fill_manual(name ="Changes in human modification",
                      values = rev(pal),
                      labels = c("high",
                                 "low")) +
    theme_minimal() +
    # xlab(glue("{expression(delta)} Human Mobility")) +
    # xlab(bquote(~Delta~"Human Mobility")) +
    # ylab(bquote(~Delta~"Space Use"))+
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          # legend.position = "none",
          plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 10, 
                                      face = "bold"),
          axis.title.x = element_text(size = 10, 
                                      face = "bold"),
          axis.ticks.x = element_line(color = "#4a4e4d"),
          text = element_text(family = "Arial", color = "#4a4e4d")) +
    labs(x = "Relative change in human mobility", y = "Standardized  \n change in area size", tag = bquote(bold("a")))
)

#- Combine plots

(comb_plot <- area_int_ce_plot + niche_int_ce_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom"))

ggsave(comb_plot, file = "figs/intra-ind_fig.png")


####----    Summarize effects Size ----####


library(emmeans)
library(bayestestR)
##-- Niche  --##

#- Get Marginal Effects at Median -#
niche_med_sg <- median(niche_out$data$sg_diff, na.rm = T)
niche_med_ghm <- median(niche_out$data$ghm_diff, na.rm = T)


# Stash df in out lists
(niche_ghm_effects <- emtrends(niche_int_mod, ~ "sg_diff", var = "ghm_diff", 
                              at = list("sg_diff" = niche_med_sg))  %>% 
  as.data.frame() %>% 
  rename("Estimate" = "ghm_diff.trend",
         "LCL" = "lower.HPD",
         "HCL" = "upper.HPD") 
)

(niche_sg_effects <- emtrends(niche_int_mod, ~ "ghm_diff", var = "sg_diff", 
                             at = list("ghm_diff" = niche_med_ghm))  %>% 
  as.data.frame() %>% 
  rename("Estimate" = "sg_diff.trend",
         "LCL" = "lower.HPD",
         "HCL" = "upper.HPD") 
)

parameters::parameters(niche_int_mod)

##-- Area  --##

#- Get Marginal Effects at Median -#
area_med_sg <- median(area_out$data$sg_diff, na.rm = T)
area_med_ghm <- median(area_out$data$ghm_diff, na.rm = T)

# Stash df in out lists
(area_ghm_effects <- emtrends(area_int_mod, ~ "sg_diff", var = "ghm_diff", 
                              at = list("sg_diff" = area_med_sg))  %>% 
  as.data.frame() %>% 
  rename("Estimate" = "ghm_diff.trend",
         "LCL" = "lower.HPD",
         "HCL" = "upper.HPD") )

(area_sg_effects <- emtrends(area_int_mod, ~ "ghm_diff", var = "sg_diff", 
                             at = list("ghm_diff" = area_med_ghm))  %>% 
  as.data.frame() %>% 
  rename("Estimate" = "sg_diff.trend",
         "LCL" = "lower.HPD",
         "HCL" = "upper.HPD") )

parameters::parameters(area_int_mod)
