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

  .wd <- '~/repositories/human_mobility_wildlife'
  .datP <- file.path(.wd,'out/single_species_models')
  .traitPF <- '/home/julietcohen/covid_movement_full_repo/raw_data/anthropause_data_sheet.csv'
  .outPF <- file.path(.wd,'out/figs/')
  
} else {

  library(docopt)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  source('src/funs/input_parse.r')
  .datP <- makePath(ag$dat)
  .outPF <- makePath(ag$out)
  .traitPF <- makePath(ag$trait)
}

#---- Initialize Environment ----#
t0 <- Sys.time()

source('src/startup.r')

suppressWarnings(
  suppressPackageStartupMessages({
    library(tidyverse)
    library(ggplot2)
    library(ggthemes)
    library(patchwork)
    library(brms)
    library(grid)
    library(emmeans)
    library(parameters)
    library(glue)
  }))

call_fun <- function(f,x,...) {
  n <- length(x)
  items <- paste(glue('x[[{1:n}]]'),collapse=',')
  eval(parse(text=glue('f({items},...)')))
}

#Source all files in the auto load funs directory
list.files('src/funs/auto',full.names=TRUE) %>%
  walk(source)

#---- Load data ----#
message("Loading trait data...")
traits <- read_csv(.traitPF) %>%
  add_row(Species = "Procyon lotor",
          class = "mammal",
          migratory = "non-migratory") %>% 
  add_row(Species = "Spilogale putorius",
          class = "mammal",
          migratory = "non-migratory") %>% 
  add_row(Species = "Sus scrofa",
          class = "mammal",
          migratory = "non-migratory") %>% 
  mutate(Species = case_when(Species == "Chen rossii" ~ "Anser rossii",
                             TRUE ~ Species))


# #-- Interaction Models --#

message('Loading interaction models...')
int_modlist <- list.files(path=file.path(.datP, "area_interactive/"), full.names = F )
int_modlist_full <- list.files( path=file.path(.datP, "area_interactive/"), full.names = T)
int_sp <- word(int_modlist, 1, sep = "_")

# #-- Additive Models --#

message('Loading additive‚ models...')
add_modlist <- list.files(path=file.path(.datP, "area_additive/"), 
                          full.names = F)
add_modlist_full <- list.files(path=file.path(.datP, "area_additive/"), 
                               full.names = T)
add_sp <- word(add_modlist, 1, sep = "_")


#check that lists are same
int_sp == add_sp

#---- Make Plots ----#

pal2 <- c("#E66100", "#5D3A9B") # 2 color pallete
pal3 <- c(pal2, "#808080") # add gray to the pallete
palnew <- c("#7552A3", "#CEBEDA")
palgray <- c("#808080", "#D3D3D3")


#-- make combined plots --#

res_out <- list()
sg_effects_out <- list()
ghm_effects_out <- list()

for(i in 1:length(int_modlist_full)){
  
  # interactive model
  if(int_modlist_full[i] != "NULL"){
    load(int_modlist_full[i]) # load model
    intmod <- out$model
    int_dat <- intmod$data %>% 
      mutate(ind = str_extract(grp, "^\\d+"))
    sp <- out$species
    
    fe <- parameters(intmod) #get fixed effects
    re <- posterior_summary(intmod, variable = c("sd_grp__Intercept", "sigma"))
    
    if(fe$pd[fe$Parameter=="b_sg_norm:ghm_scale"] < 0.95 & fe$pd[fe$Parameter=="b_sg_norm:ghm_scale"] > 0.05){ # if the interacton effect is non-sig...
      
      #... load the additive model instead.
      if(add_modlist_full[i] != "NULL"){
        load(add_modlist_full[i]) # load model
        addmod <- out$model
        add_dat <- addmod$data %>% 
          mutate(ind = str_extract(grp, "^\\d+"))
        sp <- out$species
        
        fe <- parameters(addmod) #get fixed effects
        re <- posterior_summary(addmod, variable = c("sd_grp__Intercept", "sigma"))
        
        sg_effects_out[[i]] <- parameters(addmod) %>%  #get fixed effects
          as.data.frame() %>% 
          mutate(species = sp,
                 sig_code = case_when(
                   pd > 0.05 & pd < 0.95 ~"ns_add",
                   TRUE ~ "sig_add"
                 ),
                 n_weeks = nrow(add_dat),
                 n_ind_yrs = length(unique(add_dat$grp)),
                 n_ind = length(unique(add_dat$ind))) %>% 
          filter(Parameter == "b_sg_norm")%>% 
          rename("Estimate" = "Median",
                 "LCL" = "CI_low",
                 "HCL" = "CI_high")
        
        ghm_effects_out[[i]] <- parameters(addmod) %>%  #get fixed effects
          as.data.frame() %>% 
          mutate(species = sp,
                 sig_code = case_when(
                   pd > 0.05 & pd < 0.95 ~"ns_add",
                   TRUE ~ "sig_add"
                 ),
                 n_weeks = nrow(add_dat),
                 n_ind_yrs = length(unique(add_dat$grp)),
                 n_ind = length(unique(add_dat$ind))) %>% 
          filter(Parameter == "b_ghm_scale")%>% 
          rename("Estimate" = "Median",
                 "LCL" = "CI_low",
                 "HCL" = "CI_high")
        
        
        coefdf <- tibble("species" = sp, 
                         "model" = "add",
                         
                         # SG EFFECTS
                         "sg_norm"=as.numeric(fe$Median[fe$Parameter == "b_sg_norm"]),
                         "sg_norm_lci"=fe$CI_low[fe$Parameter == "b_sg_norm"],
                         "sg_norm_uci"=fe$CI_high[fe$Parameter == "b_sg_norm"],
                         "sg_norm_pd"=fe$pd[fe$Parameter == "b_sg_norm"],
                         
                         # GHM EFFECTS
                         "ghm_scale"=as.numeric(fe$Median[fe$Parameter == "b_ghm_scale"]),
                         "ghm_scale_lci"=fe$CI_low[fe$Parameter == "b_ghm_scale"],
                         "ghm_scale_uci"=fe$CI_high[fe$Parameter == "b_ghm_scale"],
                         "ghm_scale_pd"=fe$pd[fe$Parameter == "b_ghm_scale"],
                         
                         # RANDOM EFFECTS
                         resid = re[2,1],
                         group = re[1,1]) %>% 
          mutate(sg_sign = case_when(sg_norm < 0 ~ "n",
                                     sg_norm >= 0 ~ "p"),
                 sg_sig = case_when((sg_norm_pd < 0.95 & sg_norm_pd > 0.05) ~ "N",
                                    TRUE ~ "Y"),
                 sg_display = case_when(sg_sig == "Y" ~ sg_norm,
                                        T ~ NA_real_),
                 ghm_sign = case_when(ghm_scale < 0 ~ "n",
                                      ghm_scale >= 0 ~ "p"),
                 ghm_sig = case_when((ghm_scale_pd < 0.95 & ghm_scale_pd > 0.05) ~ "N",
                                     TRUE ~ "Y"),
                 ghm_display = case_when(ghm_sig == "Y" ~ ghm_scale,
                                         T ~ NA_real_),
                 ICC = group/(group + resid),
                 var_ratio = resid/group)
        
        res_out[[i]] <- coefdf
      } # if add model is not NULL
    } else { #if the interaction CIs DONT overlap 0...
      #...then collect model output
      
      #- Standardized Effects Size -#
      
      # get median conditions
      med_sg <- median(out$data$sg_norm, na.rm = T)
      med_ghm <- median(out$data$ghm_scale, na.rm = T)
      
      #get interaction quantiled
      ghmq <- quantile(out$data$ghm_scale, probs = c(0.05, 0.95), na.rm = T)
      sgq <- quantile(out$data$sg_norm, probs = c(0.05, 0.95), na.rm = T)
      
      # Stash df in out lists
      ghm_effects_out[[i]] <- emtrends(intmod, ~ "sg_norm", var = "ghm_scale", 
                                       at = list("sg_norm" = sgq))  %>% 
        as.data.frame() %>% 
        mutate(species = sp,
               ghm_cond = c("Low", "High"),
               sig_code = case_when(
                 ghm_cond == "Low" ~"low_int",
                 ghm_cond == "High" ~ "high_int"
               ),
               n_weeks = nrow(int_dat),
               n_ind_yrs = length(unique(int_dat$grp)),
               n_ind = length(unique(int_dat$ind))) %>% 
        rename("Estimate" = "ghm_scale.trend",
               "LCL" = "lower.HPD",
               "HCL" = "upper.HPD")
      
      sg_effects_out[[i]] <- emtrends(intmod, ~ "ghm_scale", var = "sg_norm", 
                                      at = list("ghm_scale" = ghmq))  %>% 
        as.data.frame() %>% 
        mutate(species = sp,
               ghm_cond = c("Low", "High"),
               sig_code = case_when(
                 ghm_cond == "Low" ~"low_int",
                 ghm_cond == "High" ~ "high_int"
               ),
               n_weeks = nrow(int_dat),
               n_ind_yrs = length(unique(int_dat$grp)),
               n_ind = length(unique(int_dat$ind))) %>% 
        rename("Estimate" = "sg_norm.trend",
               "LCL" = "lower.HPD",
               "HCL" = "upper.HPD")
      
      
      coefdf <- tibble("species" = sp, 
                       "model" = "int",
                       
                       # SG EFFECTS
                       "sg_norm"=as.numeric(fe$Median[fe$Parameter == "b_sg_norm"]),
                       "sg_norm_lci"=fe$CI_low[fe$Parameter == "b_sg_norm"],
                       "sg_norm_uci"=fe$CI_high[fe$Parameter == "b_sg_norm"],
                       "sg_norm_pd"=fe$pd[fe$Parameter == "b_sg_norm"],
                       
                       # GHM EFFECTS
                       "ghm_scale"=as.numeric(fe$Median[fe$Parameter == "b_ghm_scale"]),
                       "ghm_scale_lci"=fe$CI_low[fe$Parameter == "b_ghm_scale"],
                       "ghm_scale_uci"=fe$CI_high[fe$Parameter == "b_ghm_scale"],
                       "ghm_scale_pd"=fe$pd[fe$Parameter == "b_ghm_scale"],
                       
                       # INTERCATION EFFECTS
                       "inter"=as.numeric(fe$Median[fe$Parameter == "b_sg_norm:ghm_scale"]),
                       "inter_lci"=fe$CI_low[fe$Parameter == "b_sg_norm:ghm_scale"],
                       "inter_uci"=fe$CI_high[fe$Parameter == "b_sg_norm:ghm_scale"],
                       "inter_pd"=fe$pd[fe$Parameter == "b_sg_norm:ghm_scale"],
                       
                       # RANDOM EFFECTS
                       resid = re[2,1],
                       group = re[1,1]) %>% 
        mutate(sg_sign = case_when(sg_norm < 0 ~ "n",
                                   sg_norm >= 0 ~ "p"),
               sg_sig = case_when((sg_norm_pd < 0.95 & sg_norm_pd > 0.05) ~ "N",
                                  TRUE ~ "Y"),
               sg_display = case_when(sg_sig == "Y" ~ sg_norm,
                                      T ~ NA_real_),
               ghm_sign = case_when(ghm_scale < 0 ~ "n",
                                    ghm_scale >= 0 ~ "p"),
               ghm_sig = case_when((ghm_scale_pd < 0.95 & ghm_scale_pd > 0.05) ~ "N",
                                   TRUE ~ "Y"),
               ghm_display = case_when(ghm_sig == "Y" ~ ghm_scale,
                                       T ~ NA_real_),
               inter_sign = case_when(inter < 0 ~ "n",
                                      inter >= 0 ~ "p"),
               inter_sig = case_when((inter_pd < 0.95 & inter_pd > 0.05) ~ "N",
                                     TRUE ~ "Y"),
               inter_display = case_when(inter_sig == "Y" ~ inter,
                                         T ~ NA_real_),
               ICC = group/(group + resid),
               var_ratio = resid/group)
      
      res_out[[i]] <- coefdf
    } # else collect the interactions
  } else {#if int is NULL...
    #... load the additive model instead.
    if(add_modlist_full[i] != "NULL"){
      load(add_modlist_full[i]) # load model
      addmod <- out$model
      add_dat <- addmod$data %>% 
        mutate(ind = str_extract(grp, "^\\d+"))
      sp <- out$species
      
      fe <- parameters(addmod) #get fixed effects
      re <- posterior_summary(addmod, variable = c("sd_grp__Intercept", "sigma"))
      
      sg_effects_out[[i]] <- parameters(addmod) %>%  #get fixed effects
        as.data.frame()  %>% 
        mutate(species = sp,
               sig_code = case_when(
                 pd > 0.05 & pd < 0.95 ~"ns_add",
                 TRUE ~ "sig_add"
               ),
               n_weeks = nrow(add_dat),
               n_ind_yrs = length(unique(add_dat$grp)),
               n_ind = length(unique(add_dat$ind))) %>% 
        filter(Parameter == "b_sg_norm") %>% 
        rename("Estimate" = "Median",
               "LCL" = "CI_low",
               "HCL" = "CI_high")
      
      ghm_effects_out[[i]] <- parameters(addmod) %>%  #get fixed effects
        as.data.frame() %>% 
        mutate(species = sp,
               sig_code = case_when(
                 pd > 0.05 & pd < 0.95 ~"ns_add",
                 TRUE ~ "sig_add"
               ),
               n_weeks = nrow(add_dat),
               n_ind_yrs = length(unique(add_dat$grp)),
               n_ind = length(unique(add_dat$ind))) %>% 
        filter(Parameter == "b_ghm_scale")%>% 
        rename("Estimate" = "Median",
               "LCL" = "CI_low",
               "HCL" = "CI_high")
      
      
      coefdf <- tibble("species" = sp, 
                       "model" = "add",
                       
                       # SG EFFECTS
                       "sg_norm"=as.numeric(fe$Median[fe$Parameter == "b_sg_norm"]),
                       "sg_norm_lci"=fe$CI_low[fe$Parameter == "b_sg_norm"],
                       "sg_norm_uci"=fe$CI_high[fe$Parameter == "b_sg_norm"],
                       "sg_norm_pd"=fe$pd[fe$Parameter == "b_sg_norm"],
                       
                       # GHM EFFECTS
                       "ghm_scale"=as.numeric(fe$Median[fe$Parameter == "b_ghm_scale"]),
                       "ghm_scale_lci"=fe$CI_low[fe$Parameter == "b_ghm_scale"],
                       "ghm_scale_uci"=fe$CI_high[fe$Parameter == "b_ghm_scale"],
                       "ghm_scale_pd"=fe$pd[fe$Parameter == "b_ghm_scale"],
                       
                       # RANDOM EFFECTS
                       resid = re[2,1],
                       group = re[1,1]) %>% 
        mutate(sg_sign = case_when(sg_norm < 0 ~ "n",
                                   sg_norm >= 0 ~ "p"),
               sg_sig = case_when((sg_norm_pd < 0.95 & sg_norm_pd > 0.05) ~ "N",
                                  TRUE ~ "Y"),
               sg_display = case_when(sg_sig == "Y" ~ sg_norm,
                                      T ~ NA_real_),
               ghm_sign = case_when(ghm_scale < 0 ~ "n",
                                    ghm_scale >= 0 ~ "p"),
               ghm_sig = case_when((ghm_scale_pd < 0.95 & ghm_scale_pd > 0.05) ~ "N",
                                   TRUE ~ "Y"),
               ghm_display = case_when(ghm_sig == "Y" ~ ghm_scale,
                                       T ~ NA_real_),
               ICC = group/(group + resid),
               var_ratio = resid/group)
      
      res_out[[i]] <- coefdf
    } #fi
  } #elese
  
}# i 

# combine dfs
res_out_df <- do.call("bind_rows", res_out)
write_csv(x = res_out_df, file = file.path(.outPF, glue("area_mod_summary_{Sys.Date()}.csv")))

sg_es_df <- do.call("bind_rows", sg_effects_out)
write_csv(x = sg_es_df, file = file.path(.outPF, glue("area_sg_effects_{Sys.Date()}.csv")))

ghm_es_df <- do.call("bind_rows", ghm_effects_out)
write_csv(x = ghm_es_df, file = file.path(.outPF, glue("area_ghm_effects_{Sys.Date()}.csv")))
