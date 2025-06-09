library(tidyverse)
library(brms)
library(ggthemes)
library(glue)
library(emmeans)
library(bayestestR)

.wd <- getwd()

add_mod_fp <- list.files(path = file.path(.wd, "out/intra_ind_models"), pattern = "^niche_intra_ind_add_mod_.*\\.rdata$", full.names = TRUE)
load(add_mod_fp)
add_mod <- out$mod

int_mod_fp <- list.files(path = file.path(.wd, "out/intra_ind_models"), pattern = "^niche_intra_ind_int_mod_.*\\.rdata$", full.names = TRUE)
load(int_mod_fp)
int_mod <- out$mod

# loo(add_mod, int_mod)
# waic(add_mod, int_mod, compare = T)
# 
# conditional_effects(int_mod)
# 
# pp_check(int_mod)
# pp_check(add_mod, type='error_scatter_avg')

fe <- fixef(add_mod) #get fixed effects
# add_re <- posterior_summary(out$model, variable = c("sd_grp__Intercept", "sigma"))
adddf <- tibble("species"=out$species, # grab estimates
                # SG EFFECTS
                "sg_diff"=as.numeric(fe["sg_diff", "Estimate"]),
                "sg_diff_lci"=fe["sg_diff", "Q2.5"],
                "sg_diff_uci"=fe["sg_diff", "Q97.5"],
                
                # GHM EFFECTS
                "ghm_diff"=as.numeric(fe["ghm_diff", "Estimate"]),
                "ghm_diff_lci"=fe["ghm_diff", "Q2.5"],
                "ghm_diff_uci"=fe["ghm_diff", "Q97.5"],
                
                # RANDOM EFFECTS
                # add_resid = add_re[2,1],
                # add_group = add_re[1,1]
) %>% 
  mutate(sg_sign = case_when(sg_diff < 0 ~ "n",
                             sg_diff >= 0 ~ "p"),
         sg_sig = case_when((sg_diff_lci < 0 & 0 < sg_diff_uci) ~ "N",
                            TRUE ~ "Y"),
         sg_display = case_when(sg_sig == "Y" ~ sg_diff,
                                T ~ NA_real_),
         code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
                                 sg_sig == "Y" & sg_sign == "p" ~ "Pos",
                                 sg_sig == "N" ~ "Non Sig"), 
                       levels=c("Neg", "Pos", "Non Sig")),
         
         ghm_sign = case_when(ghm_diff < 0 ~ "n",
                              ghm_diff >= 0 ~ "p"),
         ghm_sig = case_when((ghm_diff_lci < 0 & 0 < ghm_diff_uci) ~ "N",
                             TRUE ~ "Y"),
         ghm_display = case_when(ghm_sig == "Y" ~ ghm_diff,
                                 T ~ NA_real_),
         code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
                                 ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
                                 ghm_sig == "N" ~ "Non Sig"), 
                       levels=c("Neg", "Pos", "Non Sig")),
         
         # add_ICC = add_group/(add_group + add_resid),
         # add_var_ratio = add_resid/add_group
  )

# adddf


pal2 <- c("#E66100", "#5D3A9B") # 2 color pallete
pal3 <- c(pal2, "#808080") # add gray to the pallete
palnew <- c("#7552A3", "#CEBEDA")
palgray <- c("#808080", "#D3D3D3")



ce_sg <- conditional_effects(x=add_mod,
                             effects = "sg_diff",
                             re_formula = NA)
(sg_ce_plot <-  plot(ce_sg, 
                     plot = F,
                     rug = F,
                     line_args = list("se" = T,
                                      "color" = "black",
                                      "fill" = "gray"))[[1]] + 
    # scale_color_manual(values = palnew[3])+         
    theme_tufte() +
    # xlab(glue("{expression(delta)} Human Mobility")) +
    xlab(bquote(~Delta~"Human Mobility")) +
    ylab(bquote(~Delta~"Niche Breadth"))+
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    theme(
      axis.line = element_line(size = .5),
          # axis.text = element_blank(),
          # axis.ticks = element_blank(),
          # axis.title = element_blank(),
          aspect.ratio = 1,
          # text = element_text(family = "Roboto", size = 20)
    ))
ggsave(filename = file.path(.wd, "out/niche_intra_ind_sg.png"), sg_ce_plot,
       width = 6, height = 6)


ce_ghm <- conditional_effects(x=add_mod,
                             effects = "ghm_diff",
                             re_formula = NA)
(ghm_ce_plot <-  plot(ce_ghm, 
                     plot = F,
                     rug = F,
                     line_args = list("se" = T,
                                      "color" = "black",
                                      "fill" = "gray"))[[1]] + 
    # scale_color_manual(values = palnew[3])+         
    theme_tufte() +
    # xlab(glue("{expression(delta)} Human Mobility")) +
    xlab(bquote(~Delta~"Human Modification")) +
    ylab(bquote(~Delta~"Niche Breadth"))+
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    theme(axis.line = element_line(size = .5),
          # axis.text = element_blank(),
          # axis.ticks = element_blank(),
          # axis.title = element_blank(),
          aspect.ratio = 1
          # text = element_text(family = "Roboto", size = 20)
    ))
ggsave(filename = file.path(.wd, "out/niche_intra_ind_ghm.png"), ghm_ce_plot,
       width = 6, height = 6)

#--- INTERACTION MODELS

# get observed quantiles of ghm to set "low" and "high" human mod
ghmq <- quantile(out$data$ghm_diff, probs = c(0.2, 0.8), na.rm = T)
ce_int <- conditional_effects(x=int_mod,
                              effects = "sg_diff:ghm_diff",
                              int_conditions = list(ghm_diff = ghmq),
                              re_formula = NA,
                              prob = 0.9,
                              method = "posterior_linpred")

(int_ce_plot <-  plot(ce_int, 
                      plot = F,
                      rug = F,
                      line_args = list("se" = T))[[1]] + 
    scale_color_manual(name ="model structure",
                       values = c("#F98177","#8895BF"),
                       labels = c("low modification",
                                  "high modification")) +
    scale_fill_manual(name ="model structure",
                      values = c("#F98177","#8895BF"),
                      labels = c("low modification",
                                 "high modification")) +
    theme_minimal() +
    # xlab(glue("{expression(delta)} Human Mobility")) +
    xlab(bquote(~Delta~"Human Mobility")) +
    ylab(bquote(~Delta~"Space Use"))+
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 10, 
                                      face = "bold"),
          axis.ticks.x = element_line(color = "#4a4e4d"),
          text = element_text(family = "Arial", color = "#4a4e4d")) +
    labs(x = "Change in human mobility", y = "Change in area size")
)
ggsave(filename = file.path(.wd, "out/niche_intra_ind_int.png"), int_ce_plot)


####---- Get Marginal Effects at Median ----####
niche_med_sg <- median(int_mod$data$sg_diff, na.rm = T)
niche_med_ghm <- median(int_mod$data$ghm_diff, na.rm = T)


# Stash df in out lists
(niche_ghm_effects <- emtrends(int_mod, ~ "sg_diff", var = "ghm_diff", 
                               at = list("sg_diff" = niche_med_sg))  %>% 
    as.data.frame() %>% 
    rename("Estimate" = "ghm_diff.trend",
           "LCL" = "lower.HPD",
           "HCL" = "upper.HPD") 
)

(niche_sg_effects <- emtrends(int_mod, ~ "ghm_diff", var = "sg_diff", 
                              at = list("ghm_diff" = niche_med_ghm))  %>% 
    as.data.frame() %>% 
    rename("Estimate" = "sg_diff.trend",
           "LCL" = "lower.HPD",
           "HCL" = "upper.HPD") 
)

parameters::parameters(int_mod)
