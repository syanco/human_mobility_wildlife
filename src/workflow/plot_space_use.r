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
  
  .datP <- file.path(.wd,'out/single_species_models/area')
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
  }))

#Source all files in the auto load funs directory
list.files('analysis/src/funs/auto',full.names=TRUE) %>%
  walk(source)

# theme_set(theme_eda)

#---- Local parameters ----#



#---- Load data ----#
message("Loading trait data...")
traits <- read_csv(.traitPF)

message('Loading models...')
modlist <- list.files( path=.datP, full.names=TRUE ) 

modout <- matrix(nrow = length(modlist), ncol = 10)

for(i in 1:length(modlist)){
  load(modlist[i])
  fe <- fixef(out$model)
  m <- matrix(c(out$species, 
                as.numeric(fe["sg_norm", "Estimate"]),
                fe["sg_norm", "Q2.5"],
                fe["sg_norm", "Q97.5"],
                fe["ghm", "Estimate"],
                fe["ghm", "Q2.5"],
                fe["ghm", "Q97.5"],
                fe["sg_norm:ghm", "Estimate"],
                fe["sg_norm:ghm", "Q2.5"],
                fe["sg_norm:ghm", "Q97.5"]),
              nrow = 1)
  
  modout[i,] <- m
}

modout_df <- as.data.frame(modout)
colnames(modout_df) <- c("species",
                         "sg_norm",
                         "sg_norm_lci",
                         "sg_norm_uci",
                         "ghm",
                         "ghm_lci",
                         "ghm_uci",
                         "inter",
                         "inter_lci",
                         "inter_uci")
cols <- c(2:10)
modout_df[,cols] = apply(modout_df[,cols], 2, function(x) as.numeric(x))
modout_df <- modout_df %>% 
  mutate(Species = gsub("_", " ", modout_df$species)) %>% 
  left_join(traits) %>% 
  #TODO: go over these def siwth group
  mutate(diet = case_when(Diet.PlantO >= 50 ~ "Herbivore",
                          Diet.Fruit >= 50 ~ "Frugivore",
                          Diet.Scav >= 50 ~ "Savenger",
                          Diet.Inv >= 50 ~ "Insectivore",
                          (Diet.Vend+Diet.Vect+Diet.Vfish) >= 50 ~ "Carnivore",
                          TRUE ~ "Omnivore"),
         ghm_sign = case_when(ghm < 0 ~ "n",
                              ghm >= 0 ~ "p"),
         sg_sign = case_when(sg_norm < 0 ~ "n",
                             sg_norm >= 0 ~ "p"),
         inter_sign = case_when(inter < 0 ~ "n",
                                inter >= 0 ~ "p"),
         ghm_sig = case_when((ghm_lci < 0 & 0 < ghm_uci) ~ "N",
                             TRUE ~ "Y"),
         sg_sig = case_when((sg_norm_lci < 0 & 0 < sg_norm_uci) ~ "N",
                            TRUE ~ "Y"),
         inter_sig = case_when((inter_lci < 0 & 0 < inter_uci) ~ "N",
                               TRUE ~ "Y"),
         ghm_display = case_when(ghm_sig == "Y" ~ ghm,
                                 T ~ NA_real_),
         sg_display = case_when(sg_sig == "Y" ~ sg_norm,
                                T ~ NA_real_),
         inter_display = case_when(inter_sig == "Y" ~ inter,
                                   T ~ NA_real_))


#---- Make Plot ----#

# FIG 2

pal <- c("#E66100", "#5D3A9B")

# (pName <- modout_df %>%
#     ggplot()+
#     geom_text(aes(y = species, x=-1, label = species), hjust = 0) +
#     # ggtitle("Species") +
#     theme_void() +
#     theme(legend.position = "none"))

(pMob <- modout_df %>% 
    ggplot()+
    geom_segment(aes(color = sg_sign, x = 0, xend = sg_display, y = Species, 
                     yend = Species),
                 size = 5, alpha = 0.5) +
    scale_color_manual(values = pal) +
    geom_text(aes(y = Species, x = 50, label = round(sg_display,2))) +
    ggtitle("Human Mobility") +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(face="plain", color="black", size=10, hjust = 0)))

(pGHM <- modout_df %>%
    ggplot()+
    geom_segment(aes(color = ghm_sign, x = 0, xend = ghm_display, y = Species, 
                     yend = Species),
                 size = 5, alpha = 0.5) +
    scale_color_manual(values = pal) +
    geom_text(aes(y = Species, x = 0, label = round(ghm_display,2))) +
    ggtitle("Human Modification")+
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)))
(pInter <- modout_df %>%
    ggplot()+
    geom_segment(aes(color = inter_sign, x = 0, xend = inter_display, y = Species, 
                     yend = Species),
                 size = 5, alpha = 0.5) +
    scale_color_manual(values = pal) +
    geom_text(aes(y = Species, x = 0, label = round(inter_display,2))) +
    ggtitle("Human Mobility X Human Modification")+
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)))

pMob|pGHM|pInter

# FIG 3
pal2 <- c(pal, "#808080")
modout_df <- modout_df %>% 
  filter(Species != "Ursus arctos",
         Species != "Cervus elaphus",
         Species != "Canis latrans") 
# diet
(dietMod <- modout_df %>% 
    mutate(code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
                                   ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
                                   ghm_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=ghm, y=Species, color = code))+
    geom_errorbarh(aes(xmin=ghm_lci, xmax=ghm_uci, y = Species, color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0,linetype = "dashed")+
    theme_tufte() +
    ggtitle("Human Modification") +
    ylab("Diet Guild")+
    xlab("")+
    facet_wrap(~diet, ncol = 1, strip.position="left", scales = "free_y")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
)
(dietMob<- modout_df %>% 
    mutate(code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
                                   sg_sig == "Y" & sg_sign == "p" ~ "Pos",
                                   sg_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=sg_norm, y=Species, color = code))+
    geom_errorbarh(aes(xmin=sg_norm_lci, xmax=sg_norm_uci, y = Species, color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0, linetype = "dashed")+
    theme_tufte() +
    ggtitle("Human Mobility") +
    ylab("")+
    xlab("")+
    facet_wrap(~diet, ncol = 1, strip.position="left", scales = "free_y")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
)

(dietInter <- modout_df %>% 
    mutate(code = factor(case_when(inter_sig == "Y" & inter_sign == "n" ~ "Neg",
                                   inter_sig == "Y" & inter_sign == "p" ~ "Pos",
                                   inter_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=inter, y=Species, color = code))+
    geom_errorbarh(aes(xmin=inter_lci, xmax=inter_uci, y = Species, color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0, linetype = "dashed")+
    theme_tufte() +
    ggtitle("Interaction") +
    ylab("")+
    xlab("")+
    facet_wrap(~diet, ncol = 1, strip.position="left", scales = "free_y")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
)

# Migration
(migMod <- modout_df %>% 
    mutate(code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
                                   ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
                                   ghm_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=ghm, y=Species, color = code))+
    geom_errorbarh(aes(xmin=ghm_lci, xmax=ghm_uci, y = Species, color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0,linetype = "dashed")+
    theme_tufte() +
    ylab("Migratory Strategy")+
    xlab("")+
    facet_wrap(~migratory, ncol = 1, strip.position="left", scales = "free_y")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
)
(migMob<- modout_df %>% 
    mutate(code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
                                   sg_sig == "Y" & sg_sign == "p" ~ "Pos",
                                   sg_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=sg_norm, y=Species, color = code))+
    geom_errorbarh(aes(xmin=sg_norm_lci, xmax=sg_norm_uci, y = Species, color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0, linetype = "dashed")+
    theme_tufte() +
    ylab("")+
    xlab("")+
    facet_wrap(~migratory, ncol = 1, strip.position="left", scales = "free_y")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
)

(migInter <- modout_df %>% 
    mutate(code = factor(case_when(inter_sig == "Y" & inter_sign == "n" ~ "Neg",
                                   inter_sig == "Y" & inter_sign == "p" ~ "Pos",
                                   inter_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=inter, y=Species, color = code))+
    geom_errorbarh(aes(xmin=inter_lci, xmax=inter_uci, y = Species, color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0, linetype = "dashed")+
    theme_tufte() +
    ylab("")+
    xlab("")+
    facet_wrap(~migratory, ncol = 1, strip.position="left", scales = "free_y")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
)

# taxonomic Group
(taxMod <- modout_df %>% 
    mutate(code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
                                   ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
                                   ghm_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=ghm, y=Species, color = code))+
    geom_errorbarh(aes(xmin=ghm_lci, xmax=ghm_uci, y = Species, color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0,linetype = "dashed")+
    theme_tufte() +
    ylab("Species Group")+
    xlab("")+
    facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
)
(taxMob<- modout_df %>% 
    mutate(code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
                                   sg_sig == "Y" & sg_sign == "p" ~ "Pos",
                                   sg_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=sg_norm, y=Species, color = code))+
    geom_errorbarh(aes(xmin=sg_norm_lci, xmax=sg_norm_uci, y = Species, color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0, linetype = "dashed")+
    theme_tufte() +
    ylab("")+
    xlab("")+
    facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
)

(taxInter <- modout_df %>% 
    mutate(code = factor(case_when(inter_sig == "Y" & inter_sign == "n" ~ "Neg",
                                   inter_sig == "Y" & inter_sign == "p" ~ "Pos",
                                   inter_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=inter, y=Species, color = code))+
    geom_errorbarh(aes(xmin=inter_lci, xmax=inter_uci, y = Species, color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0, linetype = "dashed")+
    theme_tufte() +
    ylab("")+
    xlab("")+
    facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "none")
)

# Body Mass
(massMod <- modout_df %>% 
    # filter(Species != "Ursus arctos") %>% 
    mutate(code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
                                   ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
                                   ghm_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=ghm, y=log(BodyMass.Value), color = code))+
    geom_errorbarh(aes(xmin=ghm_lci, xmax=ghm_uci, y = log(BodyMass.Value), color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0)+
    theme_tufte() +
    ylab("log(body mass)")+
    xlab("")+
    # facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
    theme(legend.position = "none")
)
(massMob <- modout_df %>% 
    # filter(Species != "Ursus arctos") %>% 
    mutate(code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
                                   sg_sig == "Y" & sg_sign == "p" ~ "Pos",
                                   sg_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=sg_norm, y=log(BodyMass.Value), color = code))+
    geom_errorbarh(aes(xmin=sg_norm_lci, xmax=sg_norm_uci, y = log(BodyMass.Value), color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0)+
    theme_tufte() +
    ylab("")+
    xlab("")+
    # facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
    theme(legend.position = "none")
)
(massInter <- modout_df %>% 
    filter(Species != "Ursus arctos") %>% 
    mutate(code = factor(case_when(inter_sig == "Y" & inter_sign == "n" ~ "Neg",
                                   inter_sig == "Y" & inter_sign == "p" ~ "Pos",
                                   inter_sig == "N" ~ "Non Sig"), levels=c("Neg", "Pos", "Non Sig"))) %>% 
    ggplot()+
    geom_point(aes(x=inter, y=log(BodyMass.Value), color = code))+
    geom_errorbarh(aes(xmin=inter_lci, xmax=inter_uci, y = log(BodyMass.Value), color = code))+
    scale_color_manual(values = pal2)+
    geom_vline(xintercept = 0)+
    theme_tufte() +
    ylab("")+
    xlab("")+
    # facet_wrap(~species_group, ncol = 1, strip.position="left", scales = "free_y")+
    theme(legend.position = "none")
)

(dietMod|dietMob|dietInter)/(migMod|migMob|migInter)/(taxMod|taxMob|taxInter)/(massMod|massMob|massInter)

#---- Save output ---#
message(glue('Saving to {.outPF}'))


message(glue('Script complete in {diffmin(t0)} minutes'))