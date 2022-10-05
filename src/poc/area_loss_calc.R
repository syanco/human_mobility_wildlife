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



# #-- Interaction Models --#
# 
message('Loading interaction models...')
int_modlist <- list.files(path=file.path(.datP, "area/"), full.names = F )
int_modlist_full <- list.files( path=file.path(.datP, "area/"), full.names = T)
int_sp <- word(int_modlist, 1, sep = "_")

# 

# #-- SG + GHM Models --#
# 
message('Loading additive models...')
add_modlist <- list.files(path=file.path(.datP, "area_additive/"), 
                          full.names = F)
add_modlist_full <- list.files(path=file.path(.datP, "area_additive/"), 
                               full.names = T)
add_sp <- word(add_modlist, 1, sep = "_")

#check that lists are same
# TODO: might want to build in checks for scenations when not true
int_sp == add_sp

#---- Make Plots ----#

pal2 <- c("#E66100", "#5D3A9B") # 2 color pallete
pal3 <- c(pal2, "#808080") # add gray to the pallete
palnew <- c("#7552A3", "#CEBEDA")
palgray <- c("#808080", "#D3D3D3")


res_out <- list()
for(i in 1:length(add_modlist_full)){
  tmp <- list()
  
  # ADDITIVE MODEL
  load(add_modlist_full[i]) # load model 
  
  # posterior_summary(out$model)
  # get observed quantiles of ghm to set "low" and "high" human mod
  ghmq <- c(min(out$data$ghm_scale, na.rm = T), quantile(out$data$ghm_scale, probs = c(0.50), na.rm = T))
  sgq <- c(min(out$data$sg_norm, na.rm = T), quantile(out$data$sg_norm, probs = c(0.50), na.rm = T))
  # ndviq <- quantile(out$data$ndvi_scale, probs = c(0.95), na.rm = T)
  # lstq <- quantile(out$data$lst_scale, probs = c(0.95), na.rm = T)
  
  # Conditional Effects Plot for interaction
  ce_int <- conditional_effects(x=out$model, 
                                effects = c("sg_norm:ghm_scale"),
                                int_conditions = list(ghm_scale = ghmq,
                                                      sg_norm = sgq),
                                re_formula = NA) 
  
  nh <- ce_int$`sg_norm:ghm_scale`[1, "estimate__"]
  h <- ce_int$`sg_norm:ghm_scale`[4, "estimate__"]
  
  # get scaling params...
  m <- mean(out$data$log_area, na.rm = T)
  sd <- sd(out$data$log_area, na.rm = T)
  
  unscaled_nh <- exp((nh * sd) + m) 
  unscaled_h <- exp((h * sd) + m) 
  
  out_df <- data.frame(species = out$species,
                       human_size = unscaled_h,
                       nonhuman_size = unscaled_nh,
                       change_prop = unscaled_h/unscaled_nh,
                       diff = unscaled_nh - unscaled_h) %>% 
    mutate(diff_sq_km = diff/1000000)
  
  res_out[[i]] <- out_df
  
}

(out <- do.call("rbind", res_out))
