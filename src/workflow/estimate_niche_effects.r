#!/usr/bin/env Rscript 
# This script makes standardized effects estimates from the fitted the area models.
# Specifically, this estimates marginal effects (slopes) for both human mobility 
# (sg) and human modification (ghm).  It also produces conditional predictions of
# space use at specified quantiles of the moderating variable.


# TODO:
#  - Update docopt

# ==== Breezy Setup ====

'
Estimate Niche Effects.

Usage:
estimate_niche_effects.r <dat> <out> 
estimate_niche_effects.r (-h | --help)


Parameters:
  dat: path to input data. 
  out: path to output directory.

Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- getwd()
  
  .datP <- file.path(.wd,'out/single_species_models')
  .outPF <- file.path(.wd,'figs/area_fig.png')
  
  # vector of probabilities foer conditional estimates
  prob_vec <- c(0.2,0.8)
  
  
} else {
  library(docopt)
  # library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  
  source('src/funs/input_parse.r')
  
  #.list <- trimws(unlist(strsplit(ag$list,',')))
  .datP <- makePath(ag$dat)
  .outPF <- makePath(ag$out)
  
  # vector of probabilities foer conditional estimates
  prob_vec <- c(0.2,0.8)
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

# call_fun <- function(f,x,...) {
#   n <- length(x)
#   items <- paste(glue('x[[{1:n}]]'),collapse=',')
#   eval(parse(text=glue('f({items},...)')))
# }

#Source all files in the auto load funs directory
list.files('src/funs/auto',full.names=TRUE) %>%
  walk(source)


#---- Load data ----#
message("***LOADING DATA***")

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

# add_sp <- append(add_sp, "NULL", after = 24)
# add_modlist_full <- append(add_modlist_full, "NULL", after = 24)


#---- Estimate Effects ----#

# Init lists to store results
pred_out <- list()
sg_effects_out <- list()
ghm_effects_out <- list()


# Loop over models
for(i in 1:length(int_modlist_full)){
  
  # Interactive Model
  if(int_modlist_full[i] != "NULL"){
    load(int_modlist_full[i]) # load model
    intmod <- out$model # store as object
    sp <- out$species # extract sp
    
    fe <- fixef(intmod) # get fixed effects
    re <- posterior_summary(intmod, variable = c("sd_grp__Intercept", "sigma")) # get random effects
    
    if(fe["sg_norm:ghm_scale", "Q2.5"] < 0 & 0 < fe["sg_norm:ghm_scale", "Q97.5"]){ # if the interaction effect overlaps 0...
      
      #... load the additive model instead.
      if(add_modlist_full[i] != "NULL"){
        
        #-- ADDITIVE --#
        
        #- Make conditional predictions -#
        load(add_modlist_full[i]) # load model
        addmod <- out$model
        sp <- out$species
        out_add <- out
        fe_add <- fixef(out_add$model) #get fixed effects
        
        # get effects significance
        add_sg_sig <- ifelse(fe_add["sg_norm", "Q2.5"] < 0 & 0 < fe_add["sg_norm", "Q97.5"], "non-sig", "sig")
        add_ghm_sig <- ifelse(fe_add["ghm_scale", "Q2.5"] < 0 & 0 < fe_add["ghm_scale", "Q97.5"], "non-sig", "sig")
        
        # get observed quantiles of ghm to set "low" and "high" human mod
        ghmq <- quantile(out$data$ghm_scale, probs = prob_vec, na.rm = T)
        sgq <- quantile(out$data$sg_norm, probs = prob_vec, na.rm = T)
        
        # Conditional Effects Plot for interaction
        ce_add <- conditional_effects(x=addmod, 
                                      effects = c("sg_norm:ghm_scale"),
                                      int_conditions = list(ghm_scale = ghmq,
                                                            sg_norm = sgq),
                                      re_formula = NA, 
                                      plot = F) 
        
        # get scaling params...
        m <- mean(out$data$log_area, na.rm = T)
        sd <- sd(out$data$log_area, na.rm = T)
        
        # stash outcomes as df
        pred_out[[i]] <- ce_add$sg_norm %>% 
          mutate(ghm_case = c("low", "high", "low", "high"),
                 sg_case = c("low", "low", "high", "high"),
                 est_unscaled = (estimate__*sd)+m,
                 est_unscaled_exp = exp(est_unscaled),
                 species = sp, 
                 model = "additive",
                 ghm_sig = add_ghm_sig,
                 sg_sig = add_sg_sig,
                 tot_sig = case_when(ghm_sig == "sig" | sg_sig == "sig" ~ "sig",
                                     TRUE ~ "non-sig"))
        
        #- Get Marginal Effects at Median -#
        sg_effects_out[[i]] <- fixef(addmod) %>%  #get fixed effects
          as.data.frame() %>% 
          rename("LCL" = "Q2.5",
                 "HCL" = "Q97.5") %>% 
          mutate(species = sp) %>% 
          rownames_to_column( var = "variable") %>% 
          filter(variable == "sg_norm")%>% 
          mutate(uncertainty = abs(HCL-LCL))
        
        ghm_effects_out[[i]] <- fixef(addmod) %>%  #get fixed effects
          as.data.frame() %>% 
          rename("LCL" = "Q2.5",
                 "HCL" = "Q97.5") %>% 
          mutate(species = sp) %>% 
          rownames_to_column( var = "variable") %>% 
          filter(variable == "ghm_scale")%>% 
          mutate(uncertainty = abs(HCL-LCL))
        
        
      } # if add model is not NULL
    } else { #if the interaction CIs DONT overlap 0...
      
      #-- INTERACTIVE --#
      
      #- Make conditional predictions -#
      load(int_modlist_full[i]) # load model
      intmod <- out$model
      sp <- out$species
      
      # get observed quantiles of ghm to set "low" and "high" human mod
      ghmq <- as.numeric(quantile(out$data$ghm_scale, probs = prob_vec, na.rm = T))
      sgq <- as.numeric(quantile(out$data$sg_norm, probs = prob_vec, na.rm = T))
      
      # Conditional Effects Plot for interaction
      ce_int <- conditional_effects(x=intmod, 
                                    effects = c("sg_norm:ghm_scale"),
                                    int_conditions = list(ghm_scale = ghmq,
                                                          sg_norm = sgq),
                                    re_formula = NA, 
                                    plot = F) 
      
      # get scaling params...
      m <- mean(out$data$log_area, na.rm = T)
      sd <- sd(out$data$log_area, na.rm = T)
      
      # stash outcomes as df
      pred_out[[i]] <- ce_int$sg_norm %>% 
        mutate(ghm_case = c("low", "high", "low", "high"),
               sg_case = c("low", "low", "high", "high"),
               est_unscaled = (estimate__*sd)+m,
               est_unscaled_exp = exp(est_unscaled),
               species = sp, 
               model = "interactive", 
               int_sig = "sig", 
               tot_sig = "sig")
      
      
      #- Get Marginal Effects at Median -#
      med_sg <- median(out$data$sg_norm, na.rm = T)
      med_ghm <- median(out$data$ghm_scale, na.rm = T)
      
      # Stash df in out lists
      ghm_effects_out[[i]] <- emtrends(intmod, ~ "sg_norm", var = "ghm_scale", 
                                       at = list("sg_norm" = med_sg))  %>% 
        as.data.frame() %>% 
        rename("Estimate" = "ghm_scale.trend",
               "LCL" = "lower.HPD",
               "HCL" = "upper.HPD") %>% 
        mutate(uncertainty = abs(HCL-LCL),
               species = sp)
      
      sg_effects_out[[i]] <- emtrends(intmod, ~ "ghm_scale", var = "sg_norm", 
                                      at = list("ghm_scale" = med_ghm))  %>% 
        as.data.frame() %>% 
        rename("Estimate" = "sg_norm.trend",
               "LCL" = "lower.HPD",
               "HCL" = "upper.HPD")%>% 
        mutate(uncertainty = abs(HCL-LCL),
               species = sp)

      
    } # else collect the interactions
  } else {#if int is NULL...
    #...then load the additive model instead
    if(add_modlist_full[i] != "NULL"){
      
      #-- ADDITIVE --#
      
      #- Make conditional predictions -#
      load(add_modlist_full[i]) # load model
      addmod <- out$model
      sp <- out$species
      out_add <- out
      fe_add <- fixef(out_add$model) #get fixed effects
      
      # get effects significance
      add_sg_sig <- ifelse(fe_add["sg_norm", "Q2.5"] < 0 & 0 < fe_add["sg_norm", "Q97.5"], "non-sig", "sig")
      add_ghm_sig <- ifelse(fe_add["ghm_scale", "Q2.5"] < 0 & 0 < fe_add["ghm_scale", "Q97.5"], "non-sig", "sig")
      
      # get observed quantiles of ghm to set "low" and "high" human mod
      ghmq <- quantile(out$data$ghm_scale, probs = prob_vec, na.rm = T)
      sgq <- quantile(out$data$sg_norm, probs = prob_vec, na.rm = T)
      
      # Conditional Effects Plot for interaction
      ce_add <- conditional_effects(x=addmod, 
                                    effects = c("sg_norm:ghm_scale"),
                                    int_conditions = list(ghm_scale = ghmq,
                                                          sg_norm = sgq),
                                    re_formula = NA, 
                                    plot = F) 
      
      # get scaling params...
      m <- mean(out$data$log_area, na.rm = T)
      sd <- sd(out$data$log_area, na.rm = T)
      
      # stash outcomes as df
      pred_out[[i]] <- ce_add$sg_norm %>% 
        mutate(ghm_case = c("low", "high", "low", "high"),
               sg_case = c("low", "low", "high", "high"),
               est_unscaled = (estimate__*sd)+m,
               est_unscaled_exp = exp(est_unscaled),
               species = sp, 
               model = "additive",
               ghm_sig = add_ghm_sig,
               sg_sig = add_sg_sig,
               tot_sig = case_when(ghm_sig == "sig" | sg_sig == "sig" ~ "sig",
                                   TRUE ~ "non-sig"))
      
      #- Get Marginal Effects at Median -#
      sg_effects_out[[i]] <- fixef(addmod) %>%  #get fixed effects
        as.data.frame() %>% 
        rename("LCL" = "Q2.5",
               "HCL" = "Q97.5") %>% 
        mutate(species = sp) %>% 
        rownames_to_column( var = "variable") %>% 
        filter(variable == "sg_norm")%>% 
        mutate(uncertainty = abs(HCL-LCL))
      
      ghm_effects_out[[i]] <- fixef(addmod) %>%  #get fixed effects
        as.data.frame() %>% 
        rename("LCL" = "Q2.5",
               "HCL" = "Q97.5") %>% 
        mutate(species = sp) %>% 
        rownames_to_column( var = "variable") %>% 
        filter(variable == "ghm_scale")%>% 
        mutate(uncertainty = abs(HCL-LCL))
      
    } #fi
  } #else
  
}# i 


#---- Write Out ----#

#-- Predicted Effects --#
pred_out_df <- do.call("bind_rows", pred_out)
write_csv(x = pred_out_df, file = file.path(.outPF, glue("niche_change_prediction_{Sys.Date()}.csv")))
save(pred_out, file = file.path(.outPF, "area_change_predictions.rdata"))

#-- Marginal Effects of SG --#
sg_es_df <- do.call("bind_rows", sg_effects_out)
write_csv(x = sg_es_df, file = file.path(.outPF, glue("niche_sg_marginal_{Sys.Date()}.csv")))

#--Marginal Effects of GHM --#
ghm_es_df <- do.call("bind_rows", ghm_effects_out)
write_csv(x = ghm_es_df, file = file.path(.outPF, glue("niche_ghm_marginal_{Sys.Date()}.csv")))

