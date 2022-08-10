#!/usr/bin/env Rscript 
# This script implements the breezy philosophy: github.com/benscarlson/breezy

# ==== Breezy setup ====

'
Plot Space Use.  Generate table-like coefficient plots for the effects of human 
modification and human mobility on animal space use.

Usage:
plot_space_use.r <dat> <out> <trait>
plot_space_use.r (-h | --help)


Parameters:
  dat: path to input data. 
  out: path to output directory.
  trait: path to trait csv

Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- getwd()
  
  .datP <- file.path(.wd,'out/single_species_models')
  .traitPF <- file.path(.wd, 'raw_data/anthropause_data_sheet.csv')
  .outPF <- file.path(.wd,'figs/area_fig.png')
  
} else {
  library(docopt)
  # library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  
  source('analysis/src/funs/input_parse.r')
  
  #.list <- trimws(unlist(strsplit(ag$list,',')))
  .datP <- makePath(ag$dat)
  .traitPF <- makePath(ag$trait)
  .outPF <- makePath(ag$out)
}

#---- Initialize Environment ----#
t0 <- Sys.time()

source('analysis/src/startup.r')

suppressWarnings(
  suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(ggthemes)
    library(patchwork)
    library(brms)
    library(grid)
  }))

#Source all files in the auto load funs directory
list.files('analysis/src/funs/auto',full.names=TRUE) %>%
  walk(source)

# theme_set(theme_eda)

#---- Local parameters ----#



#---- Load data ----#
message("Loading trait data...")
traits <- read_csv(.traitPF)

# Interaction Models
message('Loading interaction models...')
int_modlist <- list.files( path=file.path(.datP, "niche/"), full.names=TRUE ) 
# 
# int_modout <- matrix(nrow = length(int_modlist), ncol = 10)
# 
# for(i in 1:length(int_modlist)){
#   load(int_modlist[i])
#   out$model
#   fe <- fixef(out$model)
#   m <- matrix(c(out$species, 
#                 as.numeric(fe["sg_norm", "Estimate"]),
#                 fe["sg_norm", "Q2.5"],
#                 fe["sg_norm", "Q97.5"],
#                 fe["ghm_scale", "Estimate"],
#                 fe["ghm_scale", "Q2.5"],
#                 fe["ghm_scale", "Q97.5"],
#                 fe["sg_norm:ghm_scale", "Estimate"],
#                 fe["sg_norm:ghm_scale", "Q2.5"],
#                 fe["sg_norm:ghm_scale", "Q97.5"]),
#               nrow = 1)
#   
#   int_modout[i,] <- m
# }
# 
# int_modout_df <- as.data.frame(int_modout)
# colnames(int_modout_df) <- c("species",
#                              "sg_norm",
#                              "sg_norm_lci",
#                              "sg_norm_uci",
#                              "ghm",
#                              "ghm_lci",
#                              "ghm_uci",
#                              "inter",
#                              "inter_lci",
#                              "inter_uci")
# cols <- c(2:10)
# int_modout_df[,cols] = apply(int_modout_df[,cols], 2, function(x) as.numeric(x))
# int_modout_df <- int_modout_df %>% 
#   mutate(Species = gsub("_", " ", int_modout_df$species)) %>% 
#   left_join(traits) %>% 
#   #TODO: go over these def siwth group
#   mutate(diet = case_when(Diet.PlantO >= 50 ~ "Herbivore",
#                           Diet.Fruit >= 50 ~ "Frugivore",
#                           Diet.Scav >= 50 ~ "Savenger",
#                           Diet.Inv >= 50 ~ "Insectivore",
#                           (Diet.Vend+Diet.Vect+Diet.Vfish) >= 50 ~ "Carnivore",
#                           TRUE ~ "Omnivore"),
#          ghm_sign = case_when(ghm < 0 ~ "n",
#                               ghm >= 0 ~ "p"),
#          sg_sign = case_when(sg_norm < 0 ~ "n",
#                              sg_norm >= 0 ~ "p"),
#          inter_sign = case_when(inter < 0 ~ "n",
#                                 inter >= 0 ~ "p"),
#          ghm_sig = case_when((ghm_lci < 0 & 0 < ghm_uci) ~ "N",
#                              TRUE ~ "Y"),
#          sg_sig = case_when((sg_norm_lci < 0 & 0 < sg_norm_uci) ~ "N",
#                             TRUE ~ "Y"),
#          inter_sig = case_when((inter_lci < 0 & 0 < inter_uci) ~ "N",
#                                TRUE ~ "Y"),
#          ghm_display = case_when(ghm_sig == "Y" ~ ghm,
#                                  T ~ NA_real_),
#          sg_display = case_when(sg_sig == "Y" ~ sg_norm,
#                                 T ~ NA_real_),
#          inter_display = case_when(inter_sig == "Y" ~ inter,
#                                    T ~ NA_real_))

# SG Models
message('Loading safegraph models...')

sg_modlist <- list.files( path=file.path(.datP, "niche_sg/"), full.names=TRUE ) 

# sg_modout <- matrix(nrow = length(sg_modlist), ncol = 4)
# 
# for(i in 1:length(sg_modlist)){
#   load(sg_modlist[i])
#   # out$model
#   fe <- fixef(out$model)
#   m <- matrix(c(out$species, 
#                 as.numeric(fe["sg_norm", "Estimate"]),
#                 fe["sg_norm", "Q2.5"],
#                 fe["sg_norm", "Q97.5"]),
#               nrow = 1)
#   
#   sg_modout[i,] <- m
# }
# 
# sg_modout_df <- as.data.frame(sg_modout)
# colnames(sg_modout_df) <- c("species",
#                             "sg_norm",
#                             "sg_norm_lci",
#                             "sg_norm_uci")
# cols <- c(2:4)
# sg_modout_df[,cols] = apply(sg_modout_df[,cols], 2, function(x) as.numeric(x))
# sg_modout_df <- sg_modout_df %>% 
#   mutate(Species = gsub("_", " ", sg_modout_df$species)) %>% 
#   left_join(traits) %>% 
#   #TODO: go over these def siwth group
#   mutate(diet = case_when(Diet.PlantO >= 50 ~ "Herbivore",
#                           Diet.Fruit >= 50 ~ "Frugivore",
#                           Diet.Scav >= 50 ~ "Savenger",
#                           Diet.Inv >= 50 ~ "Insectivore",
#                           (Diet.Vend+Diet.Vect+Diet.Vfish) >= 50 ~ "Carnivore",
#                           TRUE ~ "Omnivore"),
#          sg_sign = case_when(sg_norm < 0 ~ "n",
#                              sg_norm >= 0 ~ "p"),
#          sg_sig = case_when((sg_norm_lci < 0 & 0 < sg_norm_uci) ~ "N",
#                             TRUE ~ "Y"),
#          sg_display = case_when(sg_sig == "Y" ~ sg_norm,
#                                 T ~ NA_real_),
#          code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
#                                  sg_sig == "Y" & sg_sign == "p" ~ "Pos",
#                                  sg_sig == "N" ~ "Non Sig"), 
#                        levels=c("Neg", "Pos", "Non Sig")))

# SG Models
message('Loading safegraph models...')

ghm_modlist <- list.files( path=file.path(.datP, "niche_ghm/"), full.names=TRUE ) 

#---- Make Plots ----#

pal2 <- c("#E66100", "#5D3A9B") # 2 color pallete
pal3 <- c(pal2, "#808080") # add gray to the pallete
palnew <- palnew <- c("#7552A3", "#CEBEDA")
palgray <- c("#808080", "#D3D3D3")


#-- make combined plots --#

# !!!!!!!!! #
# TODO:  this just manually matches the speceis lists... 
# !!!!!!!!! #

# sg_modlist <- sg_modlist[-22]

# !!!!!!!!! !

row <- list()
pl <- c()
for(i in 1:length(sg_modlist)){
  tmp <- list()
  
  # SG ONLY
  load(sg_modlist[i]) # load model
  fe <- fixef(out$model) #get fixed effects
  sgdf <- tibble("species"=out$species, # grab estimates
                 "sg_norm"=as.numeric(fe["sg_norm", "Estimate"]),
                 "sg_norm_lci"=fe["sg_norm", "Q2.5"],
                 "sg_norm_uci"=fe["sg_norm", "Q97.5"]) %>% 
    mutate(sg_sign = case_when(sg_norm < 0 ~ "n",
                               sg_norm >= 0 ~ "p"),
           sg_sig = case_when((sg_norm_lci < 0 & 0 < sg_norm_uci) ~ "N",
                              TRUE ~ "Y"),
           sg_display = case_when(sg_sig == "Y" ~ sg_norm,
                                  T ~ NA_real_),
           code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
                                   sg_sig == "Y" & sg_sign == "p" ~ "Pos",
                                   sg_sig == "N" ~ "Non Sig"), 
                         levels=c("Neg", "Pos", "Non Sig")))
  
  # Get conditional effects
  ce_sg <- conditional_effects(x=out$model, 
                               effects = "sg_norm",
                               re_formula = NA)
  # if(sgdf$sg_sig == "N"){next}
  if(sgdf$sg_sig == "Y"){ # If the effect is significant...
    (sg_ce_plot <-  plot(ce_sg, plot = F,
                         line_args = list("se" = F,
                                          "color" = pal2[2]))[[1]] + 
       # scale_color_manual(values = palnew[3])+         
       theme_tufte() +
       xlab("Human Mobility") +
       ylab("Space Use")+
       theme(axis.line = element_line(size = .5),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             # aspect.ratio = 1
             ))
    plsg <- 1
  }else{ # ...not significant:
    (sg_ce_plot <-  plot(ce_sg, plot = F,
                         line_args = list("se" = F,
                                          "color" = "#808080"))[[1]] + 
       # scale_color_manual(values = palnew[3])+         
       theme_tufte() +
       xlab("Human Mobility") +
       ylab("Space Use")+
       theme(axis.line = element_line(size = .5),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             # aspect.ratio = 1
             ))
    plsg <- 0
  }
  
  
  # GHM ONLY
  load(ghm_modlist[i]) # load model
  fe <- fixef(out$model) #get fixed effects
  ghmdf <- tibble("species"=out$species, # grab estimates
                  "ghm_scale"=as.numeric(fe["ghm_scale", "Estimate"]),
                  "ghm_scale_lci"=fe["ghm_scale", "Q2.5"],
                  "ghm_scale_uci"=fe["ghm_scale", "Q97.5"]) %>% 
    mutate(ghm_sign = case_when(ghm_scale < 0 ~ "n",
                                ghm_scale >= 0 ~ "p"),
           ghm_sig = case_when((ghm_scale_lci < 0 & 0 < ghm_scale_uci) ~ "N",
                               TRUE ~ "Y"),
           ghm_display = case_when(ghm_sig == "Y" ~ ghm_scale,
                                   T ~ NA_real_),
           code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
                                   ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
                                   ghm_sig == "N" ~ "Non Sig"), 
                         levels=c("Neg", "Pos", "Non Sig")))
  
  # Get conditional effects
  ce_ghm <- conditional_effects(x=out$model, 
                                effects = "ghm_scale",
                                re_formula = NA)
  # if(sgdf$sg_sig == "N"){next}
  if(ghmdf$ghm_sig == "Y"){ # If the effect is significant...
    (ghm_ce_plot <-  plot(ce_ghm, plot = F,
                          line_args = list("se" = F,
                                           "color" = pal2[2]))[[1]] + 
       # scale_color_manual(values = palnew[3])+         
       theme_tufte() +
       xlab("Human Modification") +
       ylab("Space Use")+
       theme(axis.line = element_line(size = .5),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             # aspect.ratio = 1
       ))
    plghm <- 1
  }else{ # ...not significant:
    (ghm_ce_plot <-  plot(ce_ghm, plot = F,
                          line_args = list("se" = F,
                                           "color" = "#808080"))[[1]] + 
       # scale_color_manual(values = palnew[3])+         
       theme_tufte() +
       xlab("Human Modification") +
       ylab("Space Use")+
       theme(axis.line = element_line(size = .5),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             # aspect.ratio = 1
       ))
    plghm <- 0
  }
  
  # INTERACTION MODEL
  load(int_modlist[i])
  fe <- fixef(out$model)
  intdf <- tibble("species" = out$species,
                  "inter" = fe["sg_norm:ghm_scale", "Estimate"],
                  "inter_lci" = fe["sg_norm:ghm_scale", "Q2.5"],
                  "inter_uci" = fe["sg_norm:ghm_scale", "Q97.5"]) %>%
    mutate(inter_sign = case_when(inter < 0 ~ "n",
                                  inter >= 0 ~ "p"),
           inter_sig = case_when((inter_lci < 0 & 0 < inter_uci) ~ "N",
                                 TRUE ~ "Y"),
           inter_display = case_when(inter_sig == "Y" ~ inter,
                                     T ~ NA_real_),
           code = factor(case_when(inter_sig == "Y" & inter_sign == "n" ~ "Neg",
                                   inter_sig == "Y" & inter_sign == "p" ~ "Pos",
                                   inter_sig == "N" ~ "Non Sig"),
                         levels=c("Neg", "Pos", "Non Sig")))
  # get observed quantiles of ghm to set "low" and "high" human mod
  ghmq <- quantile(out$data$ghm_scale, probs = c(0.10, 0.90), na.rm = T)

  # Conditional Effects Plot for interaction
  ce_int <- conditional_effects(x=out$model,
                                effects = "sg_norm:ghm_scale",
                                int_conditions = list(ghm_scale = ghmq),
                                re_formula = NA)

  if(intdf$inter_sig == "Y"){ # If the effect is significant...
    (int_ce_plot <-  plot(ce_int, plot = FALSE,
                          line_args = list("se"=F,
                                           "size" = 4))[[1]] +
       theme_tufte() +
       scale_color_manual(values = palnew, name = "Human \n Modification",
                          labels = c("High", "Low")) +
       scale_fill_manual(values = palnew, name = "Human \n Modification",
                         labels = c("High", "Low")) +
       theme(axis.line = element_line(size = 4),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             aspect.ratio = 1,
             legend.position = "none"
             
             ))
    ggsave(filename = glue("out/niche_plots/{out$species}_niche_int.png"), int_ce_plot)
    plint <- 1
  }else{ # ...if the effect is not significant
    (int_ce_plot <-  plot(ce_int, plot = FALSE,
                          line_args = list("se"=F, "size" = 4))[[1]] +
       theme_tufte() +
       # xlab("Human Mobility") +
       # ylab("")+
       scale_color_grey(name = "Human \n Modification",
                        labels = c("High", "Low"))+
       scale_fill_grey(name = "Human \n Modification",
                       labels = c("High", "Low"))+
       theme(axis.line = element_line(size = 4),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             aspect.ratio = 1,
             legend.position = "none"))

    plint <- 0
  }
  
  
  tmp[[1]] <- wrap_elements(textGrob(glue("{out$species}"),
                                     gp = gpar(fontsize = 6)))
  tmp[[2]] <- sg_ce_plot 
  tmp[[3]] <- ghm_ce_plot
  tmp[[4]] <- int_ce_plot
  
  # +
  # plot_layout(widths=c(1,1,1), heights = c(1))
  row[[i]] <- tmp 
  # wrap_plots(tmp, widths=c(1,1,1), heights = c(1))
  if(plsg == 1 | plghm == 1 | plint ==1){pl[i] <- 1}else(pl[i] <- 0)
  
}

row_reduced <- row[which(pl==1)]
fin_plots <- do.call(c, row_reduced)

(comb <- wrap_plots(fin_plots, byrow = T, ncol = 4, guides = "collect", 
                    widths=c(3,5,5,5), heights = c(5)))
ggsave(filename='out/niche_plots/comb_08092022.pdf', width = 7, height = 11.5, 
       plot = comb)



#---- Trait Model ----#

# Load trait model
load("out/single_species_models/area_trait/size_trait_mod_2022-07-30.rdata")

conds <- make_conditions(out$model, vars = "diet")

ghmq_trait <- quantile(out$data$ghm_scale, probs = c(0.05, 0.95), na.rm = T)

sgq_trait <- quantile(out$data$sg_norm, probs = c(0.10, 0.90), na.rm = T)

ce_trait <- conditional_effects(x=out$model, 
                             effects = "sg_scale:ghm_scale",
                             int_conditions = list(ghm_scale = ghmq_trait
                                                   # sg_norm = c(-.2225,-0.2138)
                                                   ),
                             conditions = conds,
                             re_formula = NA)

plot(ce_trait,
     facet_args = list("scales" = "free_y"))[[1]]+
  # xlim(c(-.2225,1))+
  theme_tufte()+
  theme(axis.line = element_line(size = .5),
        # axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

ce_trait %>% 
  mutate(diet_sig = case_when(diet = ))

#---- Save output ---#
# message(glue('Saving to {.outPF}'))

message(glue('Script complete in {diffmin(t0)} minutes'))


#---- SCRATCH ----#




# ##-- SG Univariate --#
# 
# # coef plot
# (sg_coef <-  ggplot(sg_modout_df)+
#     geom_point(aes(x=sg_norm, y=Species, color = code))+
#     geom_errorbarh(aes(xmin=sg_norm_lci, xmax=sg_norm_uci, y = Species, 
#                        color = code))+
#     scale_color_manual(values = pal3, drop = F)+ #drop = F keeps "pos" in levels
#     geom_vline(xintercept = 0,linetype = "dashed")+
#     theme_tufte() +
#     ggtitle("Human Mobility") +
#     ylab("")+
#     xlab("")+
#     theme(legend.position = "none",
#           axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           strip.background = element_blank(),
#           strip.text.x = element_blank())+
#     facet_wrap(~Species, ncol = 1))
# 
# # conditional effects plot
# 
# sg_ce_plots <- list()
# 
# for(i in 1:length(sg_modlist)){
#   
#   #load the model, species, and data
#   load(sg_modlist[i]) 
#   
#   # # get observed quantiles of ghm to set "low" and "high" human mod
#   # ghmq <- quantile(out$data$ghm_scale, probs = c(0.10, 0.90), na.rm = T)
#   # sgq <- quantile(out$data$sg_norm, probs = c(0.25, 0.75), na.rm = T)
#   # # sg_seq <- seq(as.numeric(sgq[1]), as.numeric(sgq[2]), length.out = 100)
#   # 
#   
#   # Conditional Effects Plot for interaction
#   ce_sg <- conditional_effects(x=out$model, 
#                                effects = "sg_norm",
#                                https://yale.zoom.us/j/6547577953  # int_conditions = list(ghm_scale = ghmq),
#                                rug = T,
#                                # points = T,
#                                # spaghetti = T,
#                                # ndraws=100,
#                                re_formula = NA)
#   (sg_ce_plots[[i]] <-  plot(ce_sg, plot = FALSE)[[1]] +
#       theme_void() +
#       geom_hline(aes(yintercept=0), linetype = "dashed")+
#       labs(title = out$species,
#            subtitle = "Conditional Effects") +
#       xlab("Human Mobility") +
#       ylab("Space Use")
#     # scale_color_manual(values = pal3, name = "Human \n Modification",
#     #                    labels = c("High", "Low")) +
#     # scale_fill_manual(values = pal3, name = "Human \n Modification",
#     #                   labels = c("High", "Low"))
#     
#     # theme(lege)
#   )# geom_point(data = out$data, aes(x = sg_norm, y = log_area, size = ghm), alpha = 0.2, inherit.aes = F)  
#   
# }
# 
# plot(sg_ce_plots[[2]])
# 
# pdf("figs/int_cond_07292022.pdf")
# for(i in 1:length(fl)){plot(ce_plot_list[[i]])}
# dev.off()
# 
# #-- Coef "table" plot --#
# 
# 
# 
# # (pName <- modout_df %>%
# #     ggplot()+
# #     geom_text(aes(y = species, x=-1, label = species), hjust = 0) +
# #     # ggtitle("Species") +
# #     theme_void() +
# #     theme(legend.position = "none"))
# 
# (pMob <- modout_df %>% 
#     ggplot()+
#     geom_segment(aes(color = sg_sign, x = 0, xend = sg_display, y = Species, 
#                      yend = Species),
#                  size = 5, alpha = 0.5) +
#     scale_color_manual(values = pal) +
#     geom_text(aes(y = Species, x = 50, label = round(sg_display,2))) +
#     ggtitle("Human Mobility") +
#     theme_void() +
#     theme(legend.position = "none",
#           plot.title = element_text(hjust = 0.5),
#           axis.text.y = element_text(face="plain", color="black", size=10, hjust = 0)))
# 
# (pGHM <- modout_df %>%
#     ggplot()+
#     geom_segment(aes(color = ghm_sign, x = 0, xend = ghm_display, y = Species, 
#                      yend = Species),
#                  size = 5, alpha = 0.5) +
#     scale_color_manual(values = pal) +
#     geom_text(aes(y = Species, x = 0, label = round(ghm_display,2))) +
#     ggtitle("Human Modification")+
#     theme_void() +
#     theme(legend.position = "none",
#           plot.title = element_text(hjust = 0.5)))
# (pInter <- modout_df %>%
#     ggplot()+
#     geom_segment(aes(color = inter_sign, x = 0, xend = inter_display, y = Species, 
#                      yend = Species),
#                  size = 5, alpha = 0.5) +
#     scale_color_manual(values = pal) +
#     geom_text(aes(y = Species, x = 0, label = round(inter_display,2))) +
#     ggtitle("Human Mobility X Human Modification")+
#     theme_void() +
#     theme(legend.position = "none",
#           plot.title = element_text(hjust = 0.5)))
# 
# pMob|pGHM|pInter
# 
# 
# 
# #-- Dot and whisker coef plot --#
# #
# pal2 <- c(pal, "#808080")
# modout_df <- modout_df %>% 
#   filter(Species != "Ursus arctos",
#          Species != "Cervus elaphus",
#          Species != "Canis latrans") 
# # diet
# (dietMod <- modout_df %>% 
#     mutate(code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
#                                    ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
#                                    ghm_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=ghm, y=Species, color = code))+
#     geom_errorbarh(aes(xmin=ghm_lci, xmax=ghm_uci, y = Species, color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0,linetype = "dashed")+
#     theme_tufte() +
#     ggtitle("Human Modification") +
#     ylab("Diet Guild")+
#     xlab("")+
#     facet_wrap(~diet, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           legend.position = "none")
# )
# (dietMob<- modout_df %>% 
#     mutate(code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
#                                    sg_sig == "Y" & sg_sign == "p" ~ "Pos",
#                                    sg_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=sg_norm, y=Species, color = code))+
#     geom_errorbarh(aes(xmin=sg_norm_lci, xmax=sg_norm_uci, y = Species, color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0, linetype = "dashed")+
#     theme_tufte() +
#     ggtitle("Human Mobility") +
#     ylab("")+
#     xlab("")+
#     facet_wrap(~diet, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           legend.position = "none")
# )
# 
# (dietInter <- modout_df %>% 
#     mutate(code = factor(case_when(inter_sig == "Y" & inter_sign == "n" ~ "Neg",
#                                    inter_sig == "Y" & inter_sign == "p" ~ "Pos",
#                                    inter_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=inter, y=Species, color = code))+
#     geom_errorbarh(aes(xmin=inter_lci, xmax=inter_uci, y = Species, color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0, linetype = "dashed")+
#     theme_tufte() +
#     ggtitle("Interaction") +
#     ylab("")+
#     xlab("")+
#     facet_wrap(~diet, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           legend.position = "none")
# )
# 
# # Migration
# (migMod <- modout_df %>% 
#     mutate(code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
#                                    ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
#                                    ghm_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=ghm, y=Species, color = code))+
#     geom_errorbarh(aes(xmin=ghm_lci, xmax=ghm_uci, y = Species, color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0,linetype = "dashed")+
#     theme_tufte() +
#     ylab("Migratory Strategy")+
#     xlab("")+
#     facet_wrap(~migratory, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           legend.position = "none")
# )
# (migMob<- modout_df %>% 
#     mutate(code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
#                                    sg_sig == "Y" & sg_sign == "p" ~ "Pos",
#                                    sg_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=sg_norm, y=Species, color = code))+
#     geom_errorbarh(aes(xmin=sg_norm_lci, xmax=sg_norm_uci, y = Species, color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0, linetype = "dashed")+
#     theme_tufte() +
#     ylab("")+
#     xlab("")+
#     facet_wrap(~migratory, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           legend.position = "none")
# )
# 
# (migInter <- modout_df %>% 
#     mutate(code = factor(case_when(inter_sig == "Y" & inter_sign == "n" ~ "Neg",
#                                    inter_sig == "Y" & inter_sign == "p" ~ "Pos",
#                                    inter_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=inter, y=Species, color = code))+
#     geom_errorbarh(aes(xmin=inter_lci, xmax=inter_uci, y = Species, color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0, linetype = "dashed")+
#     theme_tufte() +
#     ylab("")+
#     xlab("")+
#     facet_wrap(~migratory, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           legend.position = "none")
# )
# 
# # taxonomic Group
# (taxMod <- modout_df %>% 
#     mutate(code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
#                                    ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
#                                    ghm_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=ghm, y=Species, color = code))+
#     geom_errorbarh(aes(xmin=ghm_lci, xmax=ghm_uci, y = Species, color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0,linetype = "dashed")+
#     theme_tufte() +
#     ylab("Species Group")+
#     xlab("")+
#     facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           legend.position = "none")
# )
# (taxMob<- modout_df %>% 
#     mutate(code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
#                                    sg_sig == "Y" & sg_sign == "p" ~ "Pos",
#                                    sg_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=sg_norm, y=Species, color = code))+
#     geom_errorbarh(aes(xmin=sg_norm_lci, xmax=sg_norm_uci, y = Species, color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0, linetype = "dashed")+
#     theme_tufte() +
#     ylab("")+
#     xlab("")+
#     facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           legend.position = "none")
# )
# 
# (taxInter <- modout_df %>% 
#     mutate(code = factor(case_when(inter_sig == "Y" & inter_sign == "n" ~ "Neg",
#                                    inter_sig == "Y" & inter_sign == "p" ~ "Pos",
#                                    inter_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=inter, y=Species, color = code))+
#     geom_errorbarh(aes(xmin=inter_lci, xmax=inter_uci, y = Species, color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0, linetype = "dashed")+
#     theme_tufte() +
#     ylab("")+
#     xlab("")+
#     facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           legend.position = "none")
# )
# 
# # Body Mass
# (massMod <- modout_df %>% 
#     # filter(Species != "Ursus arctos") %>% 
#     mutate(code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
#                                    ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
#                                    ghm_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=ghm, y=log(BodyMass.Value), color = code))+
#     geom_errorbarh(aes(xmin=ghm_lci, xmax=ghm_uci, y = log(BodyMass.Value), color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0)+
#     theme_tufte() +
#     ylab("log(body mass)")+
#     xlab("")+
#     # facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(legend.position = "none")
# )
# (massMob <- modout_df %>% 
#     # filter(Species != "Ursus arctos") %>% 
#     mutate(code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
#                                    sg_sig == "Y" & sg_sign == "p" ~ "Pos",
#                                    sg_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=sg_norm, y=log(BodyMass.Value), color = code))+
#     geom_errorbarh(aes(xmin=sg_norm_lci, xmax=sg_norm_uci, y = log(BodyMass.Value), color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0)+
#     theme_tufte() +
#     ylab("")+
#     xlab("")+
#     # facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(legend.position = "none")
# )
# (massInter <- modout_df %>% 
#     filter(Species != "Ursus arctos") %>% 
#     mutate(code = factor(case_when(inter_sig == "Y" & inter_sign == "n" ~ "Neg",
#                                    inter_sig == "Y" & inter_sign == "p" ~ "Pos",
#                                    inter_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
#     ggplot()+
#     geom_point(aes(x=inter, y=log(BodyMass.Value), color = code))+
#     geom_errorbarh(aes(xmin=inter_lci, xmax=inter_uci, y = log(BodyMass.Value), color = code))+
#     scale_color_manual(values = pal2)+
#     geom_vline(xintercept = 0)+
#     theme_tufte() +
#     ylab("")+
#     xlab("")+
#     # facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
#     theme(legend.position = "none")
# )
# 
# (dietMod|dietMob|dietInter)/(migMod|migMob|migInter)/(taxMod|taxMob|taxInter)/(massMod|massMob|massInter)
# 
# 
# #-- Conditional effects plots --#
# fl <- list.files("out/single_species_models/area/")
# ce_plot_list <- list()
# 
# pal3 <- c("#D55E00", "#009E73")
# 
# for(i in 1:length(fl)){
#   
#   #load the model, species, and data
#   load(glue("out/single_species_models/area/{fl[i]}")) 
#   
#   # get observed quantiles of ghm to set "low" and "high" human mod
#   ghmq <- quantile(out$data$ghm_scale, probs = c(0.10, 0.90), na.rm = T)
#   sgq <- quantile(out$data$sg_norm, probs = c(0.25, 0.75), na.rm = T)
#   # sg_seq <- seq(as.numeric(sgq[1]), as.numeric(sgq[2]), length.out = 100)
#   
#   
#   
#   
#   # Conditional Effects Plot for interaction
#   ce_int <- conditional_effects(x=out$model, 
#                                 effects = "sg_norm:ghm_scale",
#                                 int_conditions = list(ghm_scale = ghmq),
#                                 rug = T,
#                                 # points = T,
#                                 # spaghetti = T,
#                                 # ndraws=100,
#                                 re_formula = NA)
#   (ce_plot_list[[i]] <-  plot(ce_int, plot = FALSE)[[1]] +
#       theme_void() +
#       labs(title = out$species,
#            subtitle = "Conditional Effects") +
#       xlab("Human Mobility") +
#       ylab("Space Use")+
#       scale_color_manual(values = pal3, name = "Human \n Modification",
#                          labels = c("High", "Low")) +
#       scale_fill_manual(values = pal3, name = "Human \n Modification",
#                         labels = c("High", "Low"))
#     
#     # theme(lege)
#   )# geom_point(data = out$data, aes(x = sg_norm, y = log_area, size = ghm), alpha = 0.2, inherit.aes = F)  
#   
# }
# 
# plot(ce_plot_list[[1]])
# 
# pdf("figs/int_cond_07292022.pdf")
# for(i in 1:length(fl)){plot(ce_plot_list[[i]])}
# dev.off()