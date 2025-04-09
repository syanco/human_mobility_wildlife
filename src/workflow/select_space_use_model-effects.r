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
  .outPF <- file.path(.wd,'out/figs')
  .traitPF <- file.path(.wd, 'raw_data/anthropause_data_sheet.csv')
  
} else {
  library(docopt)
  # library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  
  source('src/funs/input_parse.r')
  
  #.list <- trimws(unlist(strsplit(ag$list,',')))
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
  }))

call_fun <- function(f,x,...) {
  n <- length(x)
  items <- paste(glue('x[[{1:n}]]'),collapse=',')
  eval(parse(text=glue('f({items},...)')))
}

#Source all files in the auto load funs directory
list.files('src/funs/auto',full.names=TRUE) %>%
  walk(source)

# theme_set(theme_eda)

#---- Local parameters ----#



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
# 
message('Loading interaction models...')
int_modlist <- list.files(path=file.path(.datP, "area_interactive/"), full.names = F )
int_modlist_full <- list.files( path=file.path(.datP, "area_interactive/"), full.names = T)
int_sp <- word(int_modlist, 1, sep = "_")

# 

# #-- Additive Models --#
# 
message('Loading additive‚ models...')
add_modlist <- list.files(path=file.path(.datP, "area_additive/"), 
                          full.names = F)
add_modlist_full <- list.files(path=file.path(.datP, "area_additive/"), 
                               full.names = T)
add_sp <- word(add_modlist, 1, sep = "_")



#check that lists are same
int_sp == add_sp

# int_sp <- append(int_sp, "NULL", after = 12)
# int_modlist_full <- append(int_modlist_full, "NULL", after = 12)
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
          as.data.frame()  
          mutate(species = sp,
                 sig_code = case_when(
                   pd > 0.05 & pd < 0.95 ~"ns_add",
                   TRUE ~ "sig_add"
                 ),
                 n_weeks = nrow(add_dat),
                 n_ind_yrs = length(unique(add_dat$grp)),
                 n_ind = length(unique(add_dat$ind))) %>% 
          filter(Parameter == "b_sg_norm")
        
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
          filter(Parameter == "b_ghm_scale")
       
        
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
      
      # #get conditional effects estimates
      # sg_ce <- conditional_effects(addmod, plot = F, effects = "sg_norm:ghm_scale",
      #                              int_conditions = list("sg_norm" = med_sg,
      #                                                    "ghm_scale" = ghmq))
      # ghm_ce <- conditional_effects(addmod, plot = F, effects = "ghm_scale:sg_norm",
      #                               int_conditions = list("ghm_scale" = med_ghm,
      #                                                     "sg_norm" = sgq))
      
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
        as.data.frame()  
      mutate(species = sp,
             sig_code = case_when(
               pd > 0.05 & pd < 0.95 ~"ns_add",
               TRUE ~ "sig_add"
             ),
             n_weeks = nrow(add_dat),
             n_ind_yrs = length(unique(add_dat$grp)),
             n_ind = length(unique(add_dat$ind))) %>% 
        filter(Parameter == "b_sg_norm")
      
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
        filter(Parameter == "b_ghm_scale")
      
      
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

#---- Standardized Effects plot ----#

pal2 <- c("#E66100", "#5D3A9B") # 2 color pallete
pal3 <- c(pal2, "#808080") # add gray to the pallete
palnew <- c("#7552A3", "#CEBEDA")
palgray <- c("#808080", "#D3D3D3")

coeff_pal <- c("ns_add" = "#808080", "sig_add" = "#E66100", "low_int" = "#CEBEDA", "high_int" ="#7552A3")


# SG

sg_es_df <- do.call("bind_rows", sg_effects_out) %>% 
  left_join(traits, by = c("species" = "Species")) %>% 
  # group by 'mammal' or 'bird' class, which is character, so need .fun = first
  mutate(species_by_class = fct_reorder(species, class, .fun = first))

(sg_coef_plot <- ggplot(sg_es_df) +
    geom_point(aes(x = Estimate, y = reorder(species, Estimate), 
                   color = sig_code), position = position_dodge2(width = 0.5)) +
    geom_linerange(aes(xmin = LCL, xmax = HCL, 
                       y = reorder(species, Estimate),
                       color = sig_code), position = position_dodge2(width = 0.5)) +
    scale_color_manual(name = "", values = coeff_pal) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    ylab("Species") +
    xlab("Standardized Conditional Effect") +
    ggtitle("Human Mobility") +
    theme_classic())



(sg_coef_plot_zoom <- ggplot(sg_es_df) +
    geom_point(aes(x = Estimate, y = reorder(species, Estimate), 
                   color = sig_code), position = position_dodge2(width = 0.5)) +
    geom_linerange(aes(xmin = LCL, xmax = HCL, 
                       y = reorder(species, Estimate),
                       color = sig_code), position = position_dodge2(width = 0.5)) +
    scale_color_manual(name = "", values = coeff_pal) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    ylab("Species") +
    xlab("Standardized Conditional Effect") +
    ggtitle("Human Mobility") +
    xlim(-1,1)+
    theme_classic())

(sg_coef_plot_facet <- ggplot(sg_es_df) +
    geom_point(aes(x = Estimate, y = reorder(species, Estimate), 
                   color = sig_code), position = position_dodge2(width = 0.5)) +
    geom_linerange(aes(xmin = LCL, xmax = HCL, 
                       y = reorder(species, Estimate),
                       color = sig_code), position = position_dodge2(width = 0.5)) +
    scale_color_manual(name = "", values = coeff_pal) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    ylab("") +
    xlab("Standardized Conditional Effect") +
    ggtitle("Human Mobility") +
    theme_classic() +
    facet_wrap(~class+species, scales = "free", ncol = 1, strip.position = "left",
               labeller = label_wrap_gen(width = 10))+
    theme(axis.text.y = element_blank(),
          # strip.background = element_blank(),
          # strip.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text = element_text(size = 5),
          NULL))

# change the defalt setting for saving PNGs to avoid the requirement for having a display
options(bitmapType = "cairo")
ggsave(sg_coef_plot_facet, file = file.path(.outPF, "facet_test.png"), width = 4, height = 12)
# GHM

ghm_es_df <- do.call("bind_rows", ghm_effects_out)

(ghm_coef_plot <- ggplot(ghm_es_df) +
    geom_point(aes(x = Estimate, y = reorder(species, Estimate), color = sig_code), 
               position = position_dodge2(width = 0.5)) +
    geom_linerange(aes(xmin = LCL, xmax = HCL, 
                       y = reorder(species, Estimate), color = sig_code),
                   position = position_dodge2(width = 0.5),) +
    scale_color_manual(name = "", values = coeff_pal) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    ylab("Species") +
    xlab("Standardized Conditional Effect") +
    ggtitle("Human Modification") +
    theme_classic())



(ghm_coef_plot_zoom <- ggplot(ghm_es_df) +
    geom_point(aes(x = Estimate, y = reorder(species, Estimate), color = sig_code), 
               position = position_dodge2(width = 0.5)) +
    geom_linerange(aes(xmin = LCL, xmax = HCL, 
                       y = reorder(species, Estimate), color = sig_code),
                   position = position_dodge2(width = 0.5),) +
    scale_color_manual(name = "", values = coeff_pal) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    ylab("Species") +
    xlab("Standardized Conditional Effect") +
    ggtitle("Human Modification") +
    xlim(c(-0.5, 0.5)) +
    theme_classic())

(coef_comb <- sg_coef_plot + ghm_coef_plot + plot_layout(guides = "collect"))
(coef_zoom_comb <- sg_coef_plot_zoom + ghm_coef_plot_zoom + plot_layout(guides = "collect"))

ggsave(coef_comb, file = file.path(.outPF, "space_use_coef.png"), width = 10, height = 6)
ggsave(coef_zoom_comb, file = file.path(.outPF, "space_use_coef_zoom.png"), width = 10, height = 6)
