#!/usr/bin/env Rscript 
# This script implements the breezy philosophy: github.com/benscarlson/breezy

# ==== Breezy setup ====

'
Plot Space Use.  Generate table-like coefficient plots for the effects of human 
modification and human mobility on animal niche breadth.

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
int_modlist <- list.files(path=file.path(.datP, "niche/"), full.names = F )
int_modlist_full <- list.files( path=file.path(.datP, "niche/"), full.names = T)
int_sp <- word(int_modlist, 1, sep = "_")

# 

# #-- SG + GHM Models --#
# 
message('Loading additive models...')
add_modlist <- list.files(path=file.path(.datP, "niche_additive/"), 
                          full.names = F)
add_modlist_full <- list.files(path=file.path(.datP, "niche_additive/"), 
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
  out_add <- out
  fe <- fixef(out_add$model) #get fixed effects
  adddf <- tibble("species"=out_add$species, # grab estimates
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
  
  
  # INTERACTION MODEL
  load(int_modlist_full[i])
  out_int <- out
  fe <- fixef(out_int$model)
  int_re <- posterior_summary(out_int$model, variable = c("sd_grp__Intercept", "sigma"))
  intdf <- tibble("species" = out_int$species, 
                  "inter" = fe["sg_norm:ghm_scale", "Estimate"],
                  "inter_lci" = fe["sg_norm:ghm_scale", "Q2.5"],
                  "inter_uci" = fe["sg_norm:ghm_scale", "Q97.5"],
                  int_resid = int_re[2,1],
                  int_group = int_re[1,1]) %>% 
    mutate(inter_sign = case_when(inter < 0 ~ "n",
                                  inter >= 0 ~ "p"),
           inter_sig = case_when((inter_lci < 0 & 0 < inter_uci) ~ "N",
                                 TRUE ~ "Y"),
           inter_display = case_when(inter_sig == "Y" ~ inter,
                                     T ~ NA_real_),
           code = factor(case_when(inter_sig == "Y" & inter_sign == "n" ~ "Neg",
                                   inter_sig == "Y" & inter_sign == "p" ~ "Pos",
                                   inter_sig == "N" ~ "Non Sig"), 
                         levels=c("Neg", "Pos", "Non Sig")),
           int_ICC = int_group/(int_group + int_resid),
           int_var_ratio = int_resid/int_group)
  
  # the following control flow will estimate the conditional change between 0 
  # humans and median humans *for that sp*.  It will prioritize the interaction 
  # model if there is a sig interaction effect, then the additive model if there 
  # is no sig interaction.  if there is no sig effect, it will not produce any 
  # results
  
  if(intdf$inter_sig == "Y"){ #if there's a sig interaction effect, generate prediction from in mod
    #get scaled values for 0 for SG and GHM
    mghm <- mean(out_int$data$ghm, na.rm = T)
    sdghm <- sd(out_int$data$ghm, na.rm = T)
    ghm_0 <- -mghm/sdghm
    
    msg <- mean((out_int$data$sg/out$data$cbg_area), na.rm = T)
    sdsg <- sd((out_int$data$sg/out$data$cbg_area), na.rm = T)
    sg_0 <- -msg/sdsg
    
    
    # posterior_summary(out$model)
    # get observed quantiles of ghm to set "low" and "high" human mod
    ghmq <- c(ghm_0, quantile(out_add$data$ghm_scale, probs = c(0.50), na.rm = T))
    sgq <- c(sg_0, quantile(out_add$data$sg_norm, probs = c(0.50), na.rm = T))
    
    # ghmq <- quantile(out_int$data$ghm_scale, probs = c(0.05, 0.50), na.rm = T)
    # sgq <- quantile(out_int$data$sg_norm, probs = c(0.05, 0.50), na.rm = T)
    # ndviq <- quantile(out$data$ndvi_scale, probs = c(0.95), na.rm = T)
    # lstq <- quantile(out$data$lst_scale, probs = c(0.95), na.rm = T)
    
    # Conditional Effects Plot for interaction
    ce_int <- conditional_effects(x=out_int$model, 
                                  effects = c("sg_norm:ghm_scale"),
                                  int_conditions = list(ghm_scale = ghmq,
                                                        sg_norm = sgq),
                                  re_formula = NA) 
    
    nh <- ce_int$`sg_norm:ghm_scale`[1, "estimate__"]
    h <- ce_int$`sg_norm:ghm_scale`[4, "estimate__"]
    h_ghm_only <- ce_int$`sg_norm:ghm_scale`[3, "estimate__"]
    h_sg_only <- ce_int$`sg_norm:ghm_scale`[2, "estimate__"]
    
    # get scaling params...
    m <- mean(out_int$data$total, na.rm = T)
    sd <- sd(out_int$data$total, na.rm = T)
    # 
    unscaled_nh <- (nh * sd) + m
    unscaled_h <- (h * sd) + m
    unscaled_h_ghm <- (h_ghm_only * sd) + m
    unscaled_h_sg <- (h_sg_only * sd) + m
    
    
    out_df <- data.frame(model = "interaction",
                         species = out_int$species,
                         human_size = exp(unscaled_h),
                         nonhuman_size = exp(unscaled_nh),
                         change_prop = exp(unscaled_h-unscaled_nh),
                         ghm_only_prop = exp(unscaled_h_ghm-unscaled_nh),
                         sg_only_prop = exp(unscaled_h_sg-unscaled_nh)) %>% 
      left_join(adddf)
    
    res_out[[i]] <- out_df
  }else{ #...else check for significant additive effects
    if(adddf$ghm_sig == "Y" | adddf$sg_sig == "Y"){ #if there's a sig additive effect proceed
      #get scaled values for 0 for SG and GHM
      mghm <- mean(out_add$data$ghm, na.rm = T)
      sdghm <- sd(out_add$data$ghm, na.rm = T)
      ghm_0 <- -mghm/sdghm
      
      msg <- mean((out_add$data$sg/out_add$data$cbg_area), na.rm = T)
      sdsg <- sd((out_add$data$sg/out_add$data$cbg_area), na.rm = T)
      sg_0 <- -msg/sdsg
      
      
      # posterior_summary(out$model)
      # get observed quantiles of ghm to set "low" and "high" human mod
      ghmq <- c(ghm_0, quantile(out_add$data$ghm_scale, probs = c(0.50), na.rm = T))
      sgq <- c(sg_0, quantile(out_add$data$sg_norm, probs = c(0.50), na.rm = T))
      # 
      # ghmq <- quantile(out_add$data$ghm_scale, probs = c(0.05, 0.50), na.rm = T)
      # sgq <- quantile(out_add$data$sg_norm, probs = c(0.05, 0.50), na.rm = T)
      
      
      # ndviq <- quantile(out$data$ndvi_scale, probs = c(0.95), na.rm = T)
      # lstq <- quantile(out$data$lst_scale, probs = c(0.95), na.rm = T)
      
      # Conditional Effects Plot for interaction
      ce_int <- conditional_effects(x=out_add$model, 
                                    effects = c("sg_norm:ghm_scale"),
                                    int_conditions = list(ghm_scale = ghmq,
                                                          sg_norm = sgq),
                                    re_formula = NA, 
                                    plot = F) 
      
      nh <- ce_int$`sg_norm:ghm_scale`[1, "estimate__"]
      h <- ce_int$`sg_norm:ghm_scale`[4, "estimate__"]
      h_ghm_only <- ce_int$`sg_norm:ghm_scale`[3, "estimate__"]
      h_sg_only <- ce_int$`sg_norm:ghm_scale`[2, "estimate__"]
      
      # get scaling params...
      m <- mean(out_add$data$total, na.rm = T)
      sd <- sd(out_add$data$total, na.rm = T)
      
      unscaled_nh <- (nh * sd) + m 
      unscaled_h <- (h * sd) + m
      unscaled_h_ghm <- (h_ghm_only * sd) + m 
      unscaled_h_sg <- (h_sg_only * sd) + m
      
      out_df <- data.frame(model = "additive",
                           species = out_add$species,
                           human_size = unscaled_h,
                           nonhuman_size = unscaled_nh,
                           change_prop = exp(unscaled_h-unscaled_nh),
                           ghm_only_prop = exp(unscaled_h_ghm-unscaled_nh),
                           sg_only_prop = exp(unscaled_h_sg-unscaled_nh)) %>% 
        left_join(adddf)
      
      res_out[[i]] <- out_df
    }
  }
}

(out <- do.call("rbind", res_out) %>% 
    arrange(change_prop))
write_csv(out, glue("out/prop_change_table_niche_{Sys.Date()}.csv"))



