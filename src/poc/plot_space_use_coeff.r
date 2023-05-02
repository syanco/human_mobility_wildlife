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
dot_sp == add_sp

append(int_sp, "NULL", after = 12)
int_modlist_full <- append(int_modlist_full, "NULL", after = 12)

mod_selected <- read_csv("out/area_mod_summary_2023-04-26.csv")

#---- Make Plots ----#

pal2 <- c("#E66100", "#5D3A9B") # 2 color pallete
pal3 <- c(pal2, "#808080") # add gray to the pallete
palnew <- c("#7552A3", "#CEBEDA")
palgray <- c("#808080", "#D3D3D3")

add_sig_pal <- c("#E66100", "#FFCFAB")

# "#E66100"
# "#ffcfab"

#-- make combined plots --#

res_out <- list()

for(i in 1:nrow(mod_selected)){
  # for(i in 1:2){
    if(mod_selected[i,]$model == "int"){ # For interactive models...
    
    load(int_modlist_full[i]) # load model
    intmod <- out$model
    sp <- out$species
    
    # fe <- fixef(intmod) #get fixed effects
    # re <- posterior_summary(intmod, variable = c("sd_grp__Intercept", "sigma"))  
    # get observed quantiles of ghm to set "low" and "high" human mod
    ghmq <- quantile(out$data$ghm_scale, probs = c(0.05, 0.95), na.rm = T)
    
    # Conditional Effects Plot for interaction
    ce_int <- conditional_effects(x=out$model, 
                                  effects = "sg_norm:ghm_scale",
                                  int_conditions = list(ghm_scale = ghmq),
                                  re_formula = NA)
    
    (int_ce_plot <-  plot(ce_int, plot = FALSE,
                          line_args = list("se"=T))[[1]] +
        theme_classic() +
        scale_color_manual(values = palnew, name = "Human \n Modification",
                           labels = c("High", "Low")) +
        scale_fill_manual(values = palnew, name = "Human \n Modification",
                          labels = c("High", "Low"))+
        xlab("Human Mobility") +
        ylab("Space Use")+
        theme(axis.line = element_line(size = .5),
              legend.position = "bottom")+
        # ggtitle(sp)+
        NULL)
    
    # Put results in list
    res_out[[length(res_out)+1]] <- wrap_elements(int_ce_plot + plot_annotation(title = sp)) + plot_layout(guide = "auto")
    
  }else{ # For additive models...
    load(add_modlist_full[i]) # load model
    addmod <- out$model
    sp <- out$species
    ce_sg <- conditional_effects(x=addmod, 
                                  effects = "sg_norm",
                                  re_formula = NA)
    
    # if both are non sig, skip
    if(mod_selected[i,]$sg_sig == "N" & mod_selected[i,]$ghm_sig == "N"){
      next
    }
    
    # SG Plot
    if(mod_selected[i,]$sg_sig == "Y"){ #if the SG effect is sig
      (sg_ce_plot <-  plot(ce_sg, plot = F,
                            line_args = list("se" = T,
                                             "color" = add_sig_pal[1],
                                             "fill" = add_sig_pal[2]))[[1]] + 
         # scale_color_manual(values = palnew[3])+         
         theme_classic() +
         xlab("Human Mobility") +
         ylab("Space Use")+
         theme(axis.line = element_line(size = .5),
               # axis.text = element_blank(),
               # axis.ticks = element_blank(),
               # axis.title = element_blank(),
               aspect.ratio = 1) +
         NULL)
    }else{ # if the SG effect is not sig
      (sg_ce_plot <-  plot(ce_sg, plot = F,
                            line_args = list("se" = T,
                                             "color" = palgray[1],
                                             "fill" = palgray[2]))[[1]] + 
         # scale_color_manual(values = palnew[3])+         
         theme_tufte() +
         xlab("Human Mobility") +
         ylab("Space Use")+
         theme(axis.line = element_line(size = .5),
               # axis.text = element_blank(),
               # axis.ticks = element_blank(),
               # axis.title = element_blank(),
               aspect.ratio = 1
         ))
    }
    
    # GHM Plot
    
    ce_ghm <- conditional_effects(x=out$model, 
                                  effects = "ghm_scale",
                                  re_formula = NA)
    
    if(mod_selected[i,]$ghm_sig == "Y"){ #if the GHM effect is sig
      (ghm_ce_plot <-  plot(ce_ghm, plot = F,
                            line_args = list("se" = T,
                                             "color" = add_sig_pal[1],
                                             "fill" = add_sig_pal[2]))[[1]] + 
         # scale_color_manual(values = palnew[3])+         
         theme_classic() +
         xlab("Human Modification") +
         ylab("Space Use")+
         theme(axis.line = element_line(size = .5),
               # axis.text = element_blank(),
               # axis.ticks = element_blank(),
               # axis.title = element_blank(),
               aspect.ratio = 1) +
         NULL)
    }else{ # if the GHM effect is not sig
      (ghm_ce_plot <-  plot(ce_ghm, plot = F,
                            line_args = list("se" = T,
                                             "color" = palgray[1],
                                             "fill" = palgray[2]))[[1]] + 
         # scale_color_manual(values = palnew[3])+         
         theme_tufte() +
         xlab("Human Modification") +
         ylab("Space Use")+
         theme(axis.line = element_line(size = .5),
               # axis.text = element_blank(),
               # axis.ticks = element_blank(),
               # axis.title = element_blank(),
               aspect.ratio = 1))
    }
    
    res_out[[length(res_out)+1]] <- wrap_elements(sg_ce_plot + ghm_ce_plot + plot_layout(widths = c(4,4), guides = "auto") + 
      plot_annotation(title = sp))
    
  } # if additive models
} # i loop



(space_conditionals <- res_out[[1]] / res_out[[2]] / res_out[[3]] / res_out[[4]] / res_out[[5]] / res_out[[6]] /
  res_out[[7]] / res_out[[8]] / res_out[[9]] / res_out[[10]] +
  #   res_out[[11]] / res_out[[12]] /
  # res_out[[13]] / res_out[[14]] / res_out[[15]] / res_out[[16]] / res_out[[7]] +
  plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "bottom")
)

ggsave(space_conditionals, file = "out/space_conditionals_full.png", width = 5, height = 35)


### SCARATCH ###
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
        
        ce_sg <- conditional_effects(x=addmod, 
                                     effects = "sg_norm",
                                     re_formula = NA)
        (sg_ce_plot <-  plot(ce_sg, plot = F,
                             line_args = list("se" = T,
                                              "color" = pal2[1],
                                              "fill" = "#ffcfab"))[[1]] + 
            # scale_color_manual(values = palnew[3])+         
            theme_tufte() +
            xlab("Human Mobility") +
            ylab("Space Use")+
            theme(axis.line = element_line(linewidth = 0.5)))
        
        ce_ghm <- conditional_effects(x=addmod, 
                                      effects = "ghm_scale",
                                      re_formula = NA)
        (ghm_ce_plot <-  plot(ce_ghm, plot = F,
                              line_args = list("se" = T,
                                               "color" = pal2[1],
                                               "fill" = "#ffcfab"))[[1]] + 
            # scale_color_manual(values = palnew[3])+         
            theme_tufte() +
            xlab("Human Modification") +
            ylab("Space Use")+
            theme(axis.line = element_line(linewidth = 0.5)))
        
        res_out[[i]] <- coefdf
      } # if add model is not NULL
    } else { #if the interaction CIs DONT overlap 0...
      #...then collect model output
      
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
write_csv(x=res_out_df, file = glue("out/area_mod_summary_{Sys.Date()}.csv"))

