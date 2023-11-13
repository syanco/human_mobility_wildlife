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
    library(emmeans)
  }))

call_fun <- function(f,x,...) {
  n <- length(x)
  items <- paste(glue('x[[{1:n}]]'),collapse=',')
  eval(parse(text=glue('f({items},...)')))
}

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
int_modlist <- list.files(path=file.path(.datP, "niche_interactive/"), full.names = F )
int_modlist_full <- list.files( path=file.path(.datP, "niche_interactive/"), full.names = T)
int_sp <- word(int_modlist, 1, sep = "_")

# 

# #-- Additive Models --#
# 
message('Loading additive models...')
add_modlist <- list.files(path=file.path(.datP, "niche_additive/"), 
                          full.names = F)
add_modlist_full <- list.files(path=file.path(.datP, "niche_additive/"), 
                               full.names = T)
add_sp <- word(add_modlist, 1, sep = "_")



#check that lists are same
int_sp == add_sp

add_sp <- append(add_sp, "NULL", after = 24)
add_modlist_full <- append(add_modlist_full, "NULL", after = 24)

#---- Make Plots ----#

pal2 <- c("#E66100", "#5D3A9B") # 2 color pallete
pal3 <- c(pal2, "#808080") # add gray to the pallete
palnew <- c("#7552A3", "#CEBEDA")
palgray <- c("#808080", "#D3D3D3")



#-- make combined plots --#

res_out <- list()
sg_effects_out <- list()
ghm_effects_out <- list()
area_effects_out <- list()

for(i in 1:length(int_modlist_full)){
  
  # interactive model
  if(int_modlist_full[i] != "NULL"){
    load(int_modlist_full[i]) # load model
    intmod <- out$model
    sp <- out$species
    
    fe <- fixef(intmod) #get fixed effects
    re <- posterior_summary(intmod, variable = c("sd_grp__Intercept", "sigma"))
    
    if(fe["sg_norm:ghm_scale", "Q2.5"] < 0 & 0 < fe["sg_norm:ghm_scale", "Q97.5"]){ # if the interacton effect overalps 0...
      
      #... load the additive model instead.
      if(add_modlist_full[i] != "NULL"){
        load(add_modlist_full[i]) # load model
        addmod <- out$model
        sp <- out$species
        
        area_effects_out[[i]] <- fixef(addmod) %>%  #get fixed effects
          as.data.frame() %>% 
          rename("LCL" = "Q2.5",
                 "HCL" = "Q97.5") %>% 
          mutate(species = sp,
                 sig_code = case_when(
                   LCL < 0 & HCL > 0 ~"ns_add",
                   TRUE ~ "sig_add"
                 )) %>% 
          rownames_to_column( var = "variable") %>% 
          filter(variable == "log_area_scale")
        
        sg_effects_out[[i]] <- fixef(addmod) %>%  #get fixed effects
          as.data.frame() %>% 
          rename("LCL" = "Q2.5",
                 "HCL" = "Q97.5") %>% 
          mutate(species = sp,
                 sig_code = case_when(
                   LCL < 0 & HCL > 0 ~"ns_add",
                   TRUE ~ "sig_add"
                 )) %>% 
          rownames_to_column( var = "variable") %>% 
          filter(variable == "sg_norm")
        
        ghm_effects_out[[i]] <- fixef(addmod) %>%  #get fixed effects
          as.data.frame() %>% 
          rename("LCL" = "Q2.5",
                 "HCL" = "Q97.5") %>% 
          mutate(species = sp,
                 sig_code = case_when(
                   LCL < 0 & HCL > 0 ~"ns_add",
                   TRUE ~ "sig_add"
                 )) %>% 
          rownames_to_column( var = "variable") %>% 
          filter(variable == "ghm_scale")
        # re <- posterior_summary(addmod, variable = c("sd_grp__Intercept", "sigma"))
        
        #- Standardized Effects Size -#
        
        # # get median conditions
        # med_sg <- median(out$data$sg_norm, na.rm = T)
        # med_ghm <- median(out$data$ghm_scale, na.rm = T)
        # 
        # #get conditional effects estimates
        # sg_ce <- conditional_effects(addmod, plot = F, effects = "sg_norm",
        #                              int_conditions = list("sg_norm" = med_sg))
        # ghm_ce <- conditional_effects(addmod, plot = F, effects = "ghm_scale",
        #                               int_conditions = list("ghm_scale" = med_ghm))
        
        
        # # Stash df in out lists
        # sg_effects_out[[i]] <- sg_ce[[1]] %>% 
        #   mutate(species = sp,
        #          sig_code = case_when(
        #            lower__ < 0 & upper__ > 0 ~"ns_add",
        #            TRUE ~ "sig_add"
        #          ))
        # ghm_effects_out[[i]] <- ghm_ce[[1]] %>% 
        #   mutate(species = sp,
        #          sig_code = case_when(
        #            lower__ < 0 & upper__ > 0 ~"ns_add",
        #            TRUE ~ "sig_add"
        #          ))
        
        #collect estimates  
        coefdf <- tibble("species"= sp, # grab estimates
                         "model" = "add",
                         
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
      
      
      area_effects_out[[i]] <- fixef(intmod) %>%  #get fixed effects
        as.data.frame() %>% 
        rename("LCL" = "Q2.5",
               "HCL" = "Q97.5") %>% 
        mutate(species = sp,
               sig_code = case_when(
                 LCL < 0 & HCL > 0 ~"ns_add",
                 TRUE ~ "sig_add"
               )) %>% 
        rownames_to_column( var = "variable") %>% 
        filter(variable == "log_area_scale")
      
      # Stash df in out lists
      ghm_effects_out[[i]] <- emtrends(intmod, ~ "sg_norm", var = "ghm_scale", 
                                       at = list("sg_norm" = sgq))  %>% 
        as.data.frame() %>% 
        mutate(species = sp,
               ghm_cond = c("Low", "High"),
               sig_code = case_when(
                 ghm_cond == "Low" ~"low_int",
                 ghm_cond == "High" ~ "high_int"
               )) %>% 
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
               ))%>% 
        rename("Estimate" = "sg_norm.trend",
               "LCL" = "lower.HPD",
               "HCL" = "upper.HPD")
      
      
      
      
      coefdf <- tibble("species" = sp, 
                       "model" = "int",
                       
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
      
      res_out[[i]] <- coefdf
    } # else collect the interactions
  } else {#if int is NULL...
    #...then load the additive model instead
    if(add_modlist_full[i] != "NULL"){
      load(add_modlist_full[i]) # load model
      addmod <- out$model
      sp <- out$species
      
      fe <- fixef(addmod) #get fixed effects
      re <- posterior_summary(addmod, variable = c("sd_grp__Intercept", "sigma"))
      
      area_effects_out[[i]] <- fixef(addmod) %>%  #get fixed effects
        as.data.frame() %>% 
        rename("LCL" = "Q2.5",
               "HCL" = "Q97.5") %>% 
        mutate(species = sp,
               sig_code = case_when(
                 LCL < 0 & HCL > 0 ~"ns_add",
                 TRUE ~ "sig_add"
               )) %>% 
        rownames_to_column( var = "variable") %>% 
        filter(variable == "log_area_scale")
      
      sg_effects_out[[i]] <- fixef(addmod) %>%  #get fixed effects
        as.data.frame() %>% 
        rename("LCL" = "Q2.5",
               "HCL" = "Q97.5") %>% 
        mutate(species = sp,
               sig_code = case_when(
                 LCL < 0 & HCL > 0 ~"ns_add",
                 TRUE ~ "sig_add"
               )) %>% 
        rownames_to_column( var = "variable") %>% 
        filter(variable == "sg_norm")
      
      ghm_effects_out[[i]] <- fixef(addmod) %>%  #get fixed effects
        as.data.frame() %>% 
        rename("LCL" = "Q2.5",
               "HCL" = "Q97.5") %>% 
        mutate(species = sp,
               sig_code = case_when(
                 LCL < 0 & HCL > 0 ~"ns_add",
                 TRUE ~ "sig_add"
               )) %>% 
        rownames_to_column( var = "variable") %>% 
        filter(variable == "ghm_scale")      
      #- Standardized Effects Size -#
      
      # # get median conditions
      # med_sg <- median(out$data$sg_norm, na.rm = T)
      # med_ghm <- median(out$data$ghm_scale, na.rm = T)
      # 
      # #get conditional effects estimates
      # sg_ce <- conditional_effects(addmod, plot = F, effects = "sg_norm",
      #                              int_conditions = list("sg_norm" = med_sg))
      # ghm_ce <- conditional_effects(addmod, plot = F, effects = "ghm_scale",
      #                               int_conditions = list("ghm_scale" = med_ghm))
      # 
      # 
      # 
      # # Stash df in out lists
      # sg_effects_out[[i]] <- sg_ce[[1]] %>% 
      #   mutate(species = sp,
      #          sig_code = case_when(
      #            lower__ < 0 & upper__ > 0 ~"ns_add",
      #            TRUE ~ "sig_add"
      #          ))
      # ghm_effects_out[[i]] <- ghm_ce[[1]] %>% 
      #   mutate(species = sp,
      #          sig_code = case_when(
      #            lower__ < 0 & upper__ > 0 ~"ns_add",
      #            TRUE ~ "sig_add"
      #          ))
      # 
      
      coefdf <- tibble("species"= sp, # grab estimates
                       "model" = "add",
                       
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
      
      res_out[[i]] <- coefdf
    } #fi
  } #elese
  
}# i 

# combine dfs
res_out_df <- do.call("bind_rows", res_out)
write_csv(x = res_out_df, file = glue("out/niche_mod_summary_{Sys.Date()}.csv"))

sg_es_df <- do.call("bind_rows", sg_effects_out) %>% 
  left_join(traits, by = c("species" = "Species")) 
write_csv(x = sg_es_df, file = glue("out/niche_sg_effects_{Sys.Date()}.csv"))

ghm_es_df <- do.call("bind_rows", ghm_effects_out)%>% 
  left_join(traits, by = c("species" = "Species"))
write_csv(x = ghm_es_df, file = glue("out/niche_ghm_effects_{Sys.Date()}.csv"))

area_es_df <- do.call("bind_rows", area_effects_out)%>% 
  left_join(traits, by = c("species" = "Species"))
write_csv(x = area_es_df, file = glue("out/niche_area_effects_{Sys.Date()}.csv"))


# #---- Standardized Effects plot ----#
# 
# pal2 <- c("#E66100", "#5D3A9B") # 2 color pallete
# pal3 <- c(pal2, "#808080") # add gray to the pallete
# palnew <- c("#7552A3", "#CEBEDA")
# palgray <- c("#808080", "#D3D3D3")
# 
# coeff_pal <- c("ns_add" = "#808080", "sig_add" = "#E66100", "low_int" = "#CEBEDA", "high_int" ="#7552A3")
# 
# 
# # SG
# 
# 
# 
# (sg_coef_plot <- ggplot(sg_es_df) +
#     geom_point(aes(x = Estimate, y = reorder(species, Estimate), 
#                    color = sig_code), position = position_dodge2(width = 0.5)) +
#     geom_linerange(aes(xmin = LCL, xmax = HCL, 
#                        y = reorder(species, Estimate),
#                        color = sig_code), position = position_dodge2(width = 0.5)) +
#     scale_color_manual(name = "", values = coeff_pal) +
#     geom_vline(aes(xintercept = 0), linetype = "dashed") +
#     ylab("Species") +
#     xlab("Standardized Conditional Effect") +
#     ggtitle("Human Mobility") +
#     theme_classic())
# 
# 
# 
# (sg_coef_plot_zoom <- ggplot(sg_es_df) +
#     geom_point(aes(x = Estimate, y = reorder(species, Estimate), 
#                    color = sig_code), position = position_dodge2(width = 0.5)) +
#     geom_linerange(aes(xmin = LCL, xmax = HCL, 
#                        y = reorder(species, Estimate),
#                        color = sig_code), position = position_dodge2(width = 0.5)) +
#     scale_color_manual(name = "", values = coeff_pal) +
#     geom_vline(aes(xintercept = 0), linetype = "dashed") +
#     ylab("Species") +
#     xlab("Standardized Conditional Effect") +
#     ggtitle("Human Mobility") +
#     xlim(-1,1)+
#     theme_classic())
# 
# (sg_coef_plot_facet <- ggplot(sg_es_df) +
#     geom_point(aes(x = Estimate, y = reorder(species, Estimate), 
#                    color = sig_code), position = position_dodge2(width = 0.5)) +
#     geom_linerange(aes(xmin = LCL, xmax = HCL, 
#                        y = reorder(species, Estimate),
#                        color = sig_code), position = position_dodge2(width = 0.5)) +
#     scale_color_manual(name = "", values = coeff_pal) +
#     geom_vline(aes(xintercept = 0), linetype = "dashed") +
#     ylab("") +
#     xlab("Standardized Conditional Effect") +
#     ggtitle("Human Mobility") +
#     theme_classic() +
#     facet_wrap(~class+species, scales = "free", ncol = 1, strip.position = "left",
#                labeller = label_wrap_gen(width = 10))+
#     theme(axis.text.y = element_blank(),
#           # strip.background = element_blank(),
#           # strip.text.x = element_blank(),
#           axis.ticks.y = element_blank(),
#           strip.text = element_text(size = 5),
#           NULL))
# ggsave(sg_coef_plot_facet, file = "out/facet_test.png", width = 4, height = 12)
# 
# 
# # GHM
# 
# (ghm_coef_plot <- ggplot(ghm_es_df) +
#     geom_point(aes(x = Estimate, y = reorder(species, Estimate), color = sig_code), 
#                position = position_dodge2(width = 0.5)) +
#     geom_linerange(aes(xmin = LCL, xmax = HCL, 
#                        y = reorder(species, Estimate), color = sig_code),
#                    position = position_dodge2(width = 0.5),) +
#     scale_color_manual(name = "", values = coeff_pal) +
#     geom_vline(aes(xintercept = 0), linetype = "dashed") +
#     ylab("Species") +
#     xlab("Standardized Conditional Effect") +
#     ggtitle("Human Modification") +
#     theme_classic())
# 
# 
# 
# (ghm_coef_plot_zoom <- ggplot(ghm_es_df) +
#     geom_point(aes(x = Estimate, y = reorder(species, Estimate), color = sig_code), 
#                position = position_dodge2(width = 0.5)) +
#     geom_linerange(aes(xmin = LCL, xmax = HCL, 
#                        y = reorder(species, Estimate), color = sig_code),
#                    position = position_dodge2(width = 0.5),) +
#     scale_color_manual(name = "", values = coeff_pal) +
#     geom_vline(aes(xintercept = 0), linetype = "dashed") +
#     ylab("Species") +
#     xlab("Standardized Conditional Effect") +
#     ggtitle("Human Modification") +
#     xlim(c(-0.5, 0.5)) +
#     theme_classic())
# 
# (coef_comb <- sg_coef_plot + ghm_coef_plot + plot_layout(guides = "collect"))
# (coef_zoom_comb <- sg_coef_plot_zoom + ghm_coef_plot_zoom + plot_layout(guides = "collect"))
# 
# ggsave(coef_comb, file = "out/space_use_coef.png", width = 10, height = 6)
# ggsave(coef_zoom_comb, file = "out/space_use_coef_zoom.png", width = 10, height = 6)

