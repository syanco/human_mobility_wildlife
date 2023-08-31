# stump script to test reversibility prediction method.  (relies on other scripts and is poorly connected)

library(tidyverse)
library(brms)
library(ggthemes)
library(glue)

conditional_effects(mod_2019)
conditional_effects(mod_2020)

pp_check(int_mod)
pp_check(int_mod, type='error_scatter_avg')


# adddf


pal2 <- c("#E66100", "#5D3A9B") # 2 color pallete
pal3 <- c(pal2, "#808080") # add gray to the pallete
# palnew <- c("#7552A3", "#CEBEDA")
palnew <- c("#8895BF", "#F98177")
palgray <- c("#808080", "#D3D3D3")

sg_quant_key <- quantile(dat$sg_norm, probs = c(0.10, 0.90), na.rm = T)

newdat <- data.frame(sg_norm = as.numeric(sg_quant_key),
                     ghm_scale = median(dat$ghm_scale, na.rm = T),
                     ndvi_scale = median(dat$ndvi_scale, na.rm = T),
                     tmax_scale = median(dat$tmax_scale, na.rm = T),
                     wk = 1, grp = "A")

p_10 <- predict(mod_2019, newdata = newdat[1,], re_formula = NA)
p_90 <- predict(mod_2019, newdata = newdat[2,], re_formula = NA)

ce_2019 <- conditional_effects(x=mod_2019,
                                      effects = "sg_norm",
                                      # conditions = data.frame(sg_norm = sg_quant_key),
                                      re_formula = NA)
ce_2020 <- conditional_effects(x=mod_2020,
                             effects = "sg_diff",
                             re_formula = NA)


# level_10 <- 
  
which(round(ce_2019$sg_norm$sg_norm, 2) == round(sg_quant_key[1], 2))
which(round(ce_2019$sg_norm$sg_norm, 1) == round(sg_quant_key[2], 1))
con_table_2019$sg_norm[101,]

(sg_ce_plot <-  plot(ce_2020, 
                     plot = F,
                     rug = F,
                     line_args = list("se" = T,
                                      "color" = "black",
                                      "fill" = "gray"))[[1]] + 
    # scale_color_manual(values = palnew[3])+         
    theme_tufte() +
    # xlab(glue("{expression(delta)} Human Mobility")) +
    xlab(bquote(~Delta~"Human Mobility")) +
    ylab(bquote("Space Use"))+
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = p_10[1,1])) +
    geom_hline(aes(yintercept = p_90[1,1])) +
    theme(axis.line = element_line(size = .5),
          # axis.text = element_blank(),
          # axis.ticks = element_blank(),
          # axis.title = element_blank(),
          aspect.ratio = 1,
          # text = element_text(family = "Roboto", size=20)
    ))
ggsave(filename = glue("out/area_intra_ind_sg.png"), sg_ce_plot,
       width = 6, height = 6)


ce_ghm <- conditional_effects(x=int_mod,
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
    ylab(bquote(~Delta~"Space Use"))+
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    theme(axis.line = element_line(size = .5),
          # axis.text = element_blank(),
          # axis.ticks = element_blank(),
          # axis.title = element_blank(),
          aspect.ratio = 1,
          # text = element_text(family = "Roboto", size=20)
    ))
ggsave(filename = glue("out/area_intra_ind_ghm.png"), ghm_ce_plot,
       width = 6, height = 6)



# get observed quantiles of ghm to set "low" and "high" human mod
ghmq <- quantile(out$data$ghm_diff, probs = c(0.10, 0.90), na.rm = T)
ce_int <- conditional_effects(x=int_mod,
                              effects = "sg_diff:ghm_diff",
                              int_conditions = list(ghm_diff = ghmq),
                              re_formula = NA)
(int_ce_plot <-  plot(ce_int, 
                      plot = F,
                      rug = F,
                      line_args = list("se" = T))[[1]] + 
    scale_color_manual(values = palnew, name = "Human \n Modification",
                       labels = c("High", "Low")) +
    scale_fill_manual(values = palnew, name = "Human \n Modification",
                      labels = c("High", "Low")) +
    theme_tufte() +
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
          axis.title = element_text(size = 11,
                                    face = "bold"),
          axis.ticks = element_line(color = "#4a4e4d"),
          axis.text = element_text(size = 12),
          text = element_text(family = "Arial", color = "#4a4e4d")) 
)
ggsave(filename = glue("out/area_intra_ind_int.png"), int_ce_plot, width = 5, height = 5)
