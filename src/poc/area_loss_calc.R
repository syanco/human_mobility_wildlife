#!/usr/bin/env Rscript 
# This script implements the breezy philosophy: github.com/benscarlson/breezy

# ==== Breezy setup ====

'
Plot Space Use.  Generate table-like coefficient plots for the effects of human 
modification and human mobility on animal space use.

Usage:
area_loss_calc.r <dat> <out> <trait>
area_loss_calc.r (-h | --help)


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

pal <- c("#ee7674", "#f4d5a4")


res_out <- list()
for(i in 1:length(add_modlist_full)){
  tmp <- list()
  
  # ADDITIVE MODEL
  load(add_modlist_full[i]) # load model 
  
  fe <- fixef(out$model) #get fixed effects
  adddf <- tibble("species"=out$species, # grab estimates
                  # SG EFFECTS
                  "sg_norm"=as.numeric(fe["sg_norm", "Estimate"]),
                  "sg_norm_lci"=fe["sg_norm", "Q2.5"],
                  "sg_norm_uci"=fe["sg_norm", "Q97.5"],
                  
                  # GHM EFFECTS
                  "ghm_scale"=as.numeric(fe["ghm_scale", "Estimate"]),
                  "ghm_scale_lci"=fe["ghm_scale", "Q2.5"],
                  "ghm_scale_uci"=fe["ghm_scale", "Q97.5"]) %>% 
    mutate(sg_sign = case_when(sg_norm < 0 ~ "n",
                               sg_norm >= 0 ~ "p"),
           sg_sig = case_when((sg_norm_lci < 0 & 0 < sg_norm_uci) ~ "N",
                              TRUE ~ "Y"),
           sg_display = case_when(sg_sig == "Y" ~ sg_norm,
                                  T ~ NA_real_),
           code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
                                   sg_sig == "Y" & sg_sign == "p" ~ "Pos",
                                   sg_sig == "N" ~ "Non Sig"), 
                         levels=c("Neg", "Pos", "Non Sig")),
           
           ghm_sign = case_when(ghm_scale < 0 ~ "n",
                                ghm_scale >= 0 ~ "p"),
           ghm_sig = case_when((ghm_scale_lci < 0 & 0 < ghm_scale_uci) ~ "N",
                               TRUE ~ "Y"),
           ghm_display = case_when(ghm_sig == "Y" ~ ghm_scale,
                                   T ~ NA_real_),
           code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
                                   ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
                                   ghm_sig == "N" ~ "Non Sig"), 
                         levels=c("Neg", "Pos", "Non Sig")))
  
  #get scaled values for 0 for SG and GHM
  mghm <- mean(out$data$ghm, na.rm = T)
  sdghm <- sd(out$data$ghm, na.rm = T)
  ghm_0 <- (0/sdghm)-mghm
  
  msg <- mean((out$data$sg/out$data$cbg_area), na.rm = T)
  sdsg <- sd((out$data$sg/out$data$cbg_area), na.rm = T)
  sg_0 <- (0/sdsg)-msg
  
  
  # posterior_summary(out$model)
  # get observed quantiles of ghm to set "low" and "high" human mod
  ghmq <- c(ghm_0, quantile(out$data$ghm_scale, probs = c(0.50), na.rm = T))
  sgq <- c(sg_0, quantile(out$data$sg_norm, probs = c(0.50), na.rm = T))
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
                       diff = unscaled_h - unscaled_nh) %>% 
    mutate(diff_sq_km = diff/1000000) %>% 
    left_join(adddf)
  
  res_out[[i]] <- out_df
  
}

(out <- do.call("rbind", res_out) %>% 
    arrange(change_prop) %>% 
    mutate(perc_change = (diff/nonhuman_size)*100)
)
write_csv(out, "out/prop_change_table_20221006.csv")


# Plot
traits <- read_csv("raw_data/anthropause_data_sheet.csv")


out %>% 
  filter(sg_sig == "Y" | ghm_sig == "Y") %>% 
  left_join(traits, by = c("species" = "Species")) %>% 
  mutate(effect = case_when(diff < 0 ~ "neg",
                            diff > 0 ~ "pos"),
         class = case_when(class == "bird" ~ "Birds",
                           class == "mammal" ~ "Mammals"),
         species = fct_reorder(species, perc_change, .desc = T)) %>% 
  ggplot()+
  geom_bar(aes(y = perc_change, x = species, fill = effect), stat = "identity", width = .5)+
  geom_hline(aes(yintercept = 0)) +
  ylab("Percent Change in Weekly Space Use") +
  coord_flip() +
  scale_fill_manual(values = pal) +
  facet_wrap(~class, scales = "free_y", nrow = 2) +
  theme_tufte()+
  theme(legend.position = "none")


