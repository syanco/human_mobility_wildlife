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
int_modlist <- list.files(path=file.path(.datP, "area_interactive/"), full.names = F )
int_modlist_full <- list.files( path=file.path(.datP, "area_interactive/"), full.names = T)
int_sp <- word(int_modlist, 1, sep = "_")

# 

# #-- Additive Models --#
# 
message('Loading additive models...')
add_modlist <- list.files(path=file.path(.datP, "area_additive/"), 
                          full.names = F)
add_modlist_full <- list.files(path=file.path(.datP, "area_additive/"), 
                               full.names = T)
add_sp <- word(add_modlist, 1, sep = "_")

# #-- Dot Models --#
# 
message('Loading intercept-only  models...')
dot_modlist <- list.files(path=file.path(.datP, "area_dot/"), 
                          full.names = F)
dot_modlist_full <- list.files(path=file.path(.datP, "area_dot/"), 
                               full.names = T)
dot_sp <- word(dot_modlist, 1, sep = "_")


#check that lists are same
# TODO: might want to build in checks for scenations when not true
int_sp == add_sp
int_sp == dot_sp

#---- Make Plots ----#

pal2 <- c("#E66100", "#5D3A9B") # 2 color pallete
pal3 <- c(pal2, "#808080") # add gray to the pallete
palnew <- c("#7552A3", "#CEBEDA")
palgray <- c("#808080", "#D3D3D3")



#-- make combined plots --#

res_out <- list()

for(i in 1:length(add_modlist_full)){
  
  # dot model
  load(dot_modlist_full[i]) # load model
  dotmod <- out$model
  sp <- out$species
  
  
  # additive model
  load(add_modlist_full[i]) # load model
  addmod <- out$model
  
  # interactive model
  load(int_modlist_full[i]) # load model
  intmod <- out$model
  
  # execute model compare
  comp <- loo(dotmod, addmod, intmod)
  
  winner <- names(which(comp$diffs[,"elpd_diff"] == 0))
  focal <- get(winner)
  
  fe <- fixef(focal) #get fixed effects
  re <- posterior_summary(focal, variable = c("sd_grp__Intercept", "sigma"))
  
  if(winner == "dotmod"){ # if dot model...
    coefdf <- tibble("species" = sp, 
                     "model" = winner,
                     "dot" = fe["Intercept", "Estimate"],
                     "dot_lci" = fe["Intercept", "Q2.5"],
                     "dot_uci" = fe["Intercept", "Q97.5"],
                     "resid" = re[2,1],
                     "group" = re[1,1]) %>% 
      mutate(dot_sign = case_when(dot < 0 ~ "n",
                                  dot >= 0 ~ "p"),
             dot_sig = case_when((dot_lci < 0 & 0 < dot_uci) ~ "N",
                                 TRUE ~ "Y"),
             dot_display = case_when(dot_sig == "Y" ~ dot,
                                     T ~ NA_real_),
             ICC = group/(group + resid),
             var_ratio = resid/group)
    
  } else {
    if(winner == "addmod"){ # if add model
      coefdf <- tibble("species"= sp, # grab estimates
                       "model" = winner,
                       
                       # SG EFFECTS
                       "sg_norm"=as.numeric(fe["sg_norm", "Estimate"]),
                       "sg_norm_lci"=fe["sg_norm", "Q2.5"],
                       "sg_norm_uci"=fe["sg_norm", "Q97.5"],
                       
                       # GHM EFFECTS
                       "ghm_scale"=as.numeric(fe["ghm_scale", "Estimate"]),
                       "ghm_scale_lci"=fe["ghm_scale", "Q2.5"],
                       "ghm_scale_uci"=fe["ghm_scale", "Q97.5"],
                       
                       # RANDOM EFFECTS
                       resid = re[2,1],
                       group = re[1,1]) %>% 
        mutate(sg_sign = case_when(sg_norm < 0 ~ "n",
                                   sg_norm >= 0 ~ "p"),
               sg_sig = case_when((sg_norm_lci < 0 & 0 < sg_norm_uci) ~ "N",
                                  TRUE ~ "Y"),
               sg_display = case_when(sg_sig == "Y" ~ sg_norm,
                                      T ~ NA_real_),
               ghm_sign = case_when(ghm_scale < 0 ~ "n",
                                    ghm_scale >= 0 ~ "p"),
               ghm_sig = case_when((ghm_scale_lci < 0 & 0 < ghm_scale_uci) ~ "N",
                                   TRUE ~ "Y"),
               ghm_display = case_when(ghm_sig == "Y" ~ ghm_scale,
                                       T ~ NA_real_),
               ICC = group/(group + resid),
               var_ratio = resid/group)
      
    } else {
      if(winner == "intmod"){ # if int model
        coefdf <- tibble("species" = sp, 
                         "model" = winner,
                         
                         # SG EFFECTS
                         "sg_norm"=as.numeric(fe["sg_norm", "Estimate"]),
                         "sg_norm_lci"=fe["sg_norm", "Q2.5"],
                         "sg_norm_uci"=fe["sg_norm", "Q97.5"],
                         
                         # GHM EFFECTS
                         "ghm_scale"=as.numeric(fe["ghm_scale", "Estimate"]),
                         "ghm_scale_lci"=fe["ghm_scale", "Q2.5"],
                         "ghm_scale_uci"=fe["ghm_scale", "Q97.5"],
                         
                         # INTERCATION EFFECTS
                         "inter" = fe["sg_norm:ghm_scale", "Estimate"],
                         "inter_lci" = fe["sg_norm:ghm_scale", "Q2.5"],
                         "inter_uci" = fe["sg_norm:ghm_scale", "Q97.5"],
                         
                         # RANDOM EFFECTS
                         resid = re[2,1],
                         group = re[1,1]) %>% 
          mutate(sg_sign = case_when(sg_norm < 0 ~ "n",
                                     sg_norm >= 0 ~ "p"),
                 sg_sig = case_when((sg_norm_lci < 0 & 0 < sg_norm_uci) ~ "N",
                                    TRUE ~ "Y"),
                 sg_display = case_when(sg_sig == "Y" ~ sg_norm,
                                        T ~ NA_real_),
                 ghm_sign = case_when(ghm_scale < 0 ~ "n",
                                      ghm_scale >= 0 ~ "p"),
                 ghm_sig = case_when((ghm_scale_lci < 0 & 0 < ghm_scale_uci) ~ "N",
                                     TRUE ~ "Y"),
                 ghm_display = case_when(ghm_sig == "Y" ~ ghm_scale,
                                         T ~ NA_real_),
                 inter_sign = case_when(inter < 0 ~ "n",
                                        inter >= 0 ~ "p"),
                 inter_sig = case_when((inter_lci < 0 & 0 < inter_uci) ~ "N",
                                       TRUE ~ "Y"),
                 inter_display = case_when(inter_sig == "Y" ~ inter,
                                           T ~ NA_real_),
                 ICC = group/(group + resid),
                 var_ratio = resid/group)
      }
    }
  }
  
  res_out[[i]] <- coefdf
}# i 


# combine dfs
res_out_df <- do.call("bind_rows", res_out)
write_csv(x=res_out_df, file = glue("out/area_mod_summary_{Sys.Date()}.csv"))

