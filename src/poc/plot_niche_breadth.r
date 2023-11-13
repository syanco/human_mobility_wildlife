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
message('Loading interaction models...')
int_modlist <- list.files(path=file.path(.datP, "niche/"), 
                          pattern = "*10-14.rdata|*10-15.rdata", full.names = F)
int_modlist_full <- list.files(path=file.path(.datP, "niche/"), 
                               pattern = "*10-14.rdata|*10-15.rdata", 
                               full.names = T)
int_sp <- word(int_modlist, 1, sep = "_")

# 

# #-- SG + GHM Models --#
# 
message('Loading additive models...')
add_modlist <- list.files(path=file.path(.datP, "niche_additive/"), 
                          pattern = "*10-14.rdata|*10-15.rdata", full.names = F)
add_modlist_full <- list.files(path=file.path(.datP, "niche_additive/"), 
                               pattern = "*10-14.rdata|*10-15.rdata", 
                               full.names = T)
add_sp <- word(add_modlist, 1, sep = "_")


message("Loading area controlled models...")
cont_modlist <- list.files(path=file.path(.datP, "niche_controlled/"), 
                           pattern = "*10-14.rdata|*10-15.rdata", 
                           full.names = F)
cont_modlist_full <- list.files(path=file.path(.datP, "niche_controlled/"), 
                                pattern = "*10-14.rdata|*10-15.rdata",
                                full.names = T)
cont_sp <- word(cont_modlist, 1, sep = "_")
#check that lists are same
# TODO: might want to build in checks for scenations when not true
int_sp == add_sp
int_sp == cont_sp



#---- Make Plots ----#

# Define some color palettes
pal2 <- c("#E66100", "#5D3A9B") # 2 color palette
pal3 <- c(pal2, "#808080") # add gray to the palette
palnew <- palnew <- c("#7552A3", "#CEBEDA")
palgray <- c("#808080", "#D3D3D3")


#-- make combined plots --#

res_out <- list()
row <- list()
pl <- c()
for(i in 1:length(add_modlist_full)){
  tmp <- list()
  
  # ADDITIVE MODELS
  load(add_modlist_full[i]) # load model
  fe <- fixef(out$model) #get fixed effects
  adddf <- tibble("species"=out$species, # grab estimates
                  
                  # SG OUT
                  "add_sg_norm"=as.numeric(fe["sg_norm", "Estimate"]),
                  "add_sg_norm_lci"=fe["sg_norm", "Q2.5"],
                  "add_sg_norm_uci"=fe["sg_norm", "Q97.5"],
                  
                  # GHM OUT
                  "add_ghm_scale"=as.numeric(fe["ghm_scale", "Estimate"]),
                  "add_ghm_scale_lci"=fe["ghm_scale", "Q2.5"],
                  "add_ghm_scale_uci"=fe["ghm_scale", "Q97.5"]) %>% 
    mutate(add_sg_sign = case_when(
      add_sg_norm < 0 ~ "n",
      add_sg_norm >= 0 ~ "p"),
      add_sg_sig = case_when(
        (add_sg_norm_lci < 0 & 0 < add_sg_norm_uci) ~ "N",
        TRUE ~ "Y"),
      add_sg_display = case_when(add_sg_sig == "Y" ~ add_sg_norm,
                                 T ~ NA_real_),
      add_sg_code = factor(case_when(
        add_sg_sig == "Y" & add_sg_sign == "n" ~ "Neg",
        add_sg_sig == "Y" & add_sg_sign == "p" ~ "Pos",
        add_sg_sig == "N" ~ "Non Sig"), 
        levels=c("Neg", "Pos", "Non Sig")),
      add_ghm_sign = case_when(add_ghm_scale < 0 ~ "n",
                               add_ghm_scale >= 0 ~ "p"),
      add_ghm_sig = case_when(
        (add_ghm_scale_lci < 0 & 0 < add_ghm_scale_uci) ~ "N",
        TRUE ~ "Y"),
      add_ghm_display = case_when(add_ghm_sig == "Y" ~ add_ghm_scale,
                                  T ~ NA_real_),
      add_ghm_code = factor(case_when(
        add_ghm_sig == "Y" & add_ghm_sign == "n" ~ "Neg",
        add_ghm_sig == "Y" & add_ghm_sign == "p" ~ "Pos",
        add_ghm_sig == "N" ~ "Non Sig"), 
        levels=c("Neg", "Pos", "Non Sig")))
  
  # CONTROLLED MODEL
  load(cont_modlist_full[i]) # load model
  fe <- fixef(out$model) #get fixed effects
  # re <- posterior_summary(out$model, variable = c("sd_grp__Intercept", "sigma"))
  contdf <- tibble("species"=out$species, # grab estimates
                   
                   # AREA EFFECTS
                   "cont_area"=as.numeric(fe["log_area_scale", "Estimate"]),
                   "cont_area_lci"=fe["log_area_scale", "Q2.5"],
                   "cont_area_uci"=fe["log_area_scale", "Q97.5"],
                   
                   # SG EFFECTS
                   "cont_sg_norm"=as.numeric(fe["sg_norm", "Estimate"]),
                   "cont_sg_norm_lci"=fe["sg_norm", "Q2.5"],
                   "cont_sg_norm_uci"=fe["sg_norm", "Q97.5"],
                   
                   # GHM EFFECTS
                   "cont_ghm_scale"=as.numeric(fe["ghm_scale", "Estimate"]),
                   "cont_ghm_scale_lci"=fe["ghm_scale", "Q2.5"],
                   "cont_ghm_scale_uci"=fe["ghm_scale", "Q97.5"]) %>% 
    mutate(cont_area_sign = case_when(cont_area < 0 ~ "n",
                                      cont_area >= 0 ~ "p"),
           cont_area_sig = case_when((cont_area_lci < 0 & 0 < cont_area_uci) ~ "N",
                                     TRUE ~ "Y"),
           cont_area_display = case_when(cont_area_sig == "Y" ~ cont_area,
                                         T ~ NA_real_),
           cont_area_code = factor(case_when(
             cont_area_sig == "Y" & cont_area_sign == "n" ~ "Neg",
             cont_area_sig == "Y" & cont_area_sign == "p" ~ "Pos",
             cont_area_sig == "N" ~ "Non Sig"), 
             levels=c("Neg", "Pos", "Non Sig")),
           
           cont_sg_sign = case_when(cont_sg_norm < 0 ~ "n",
                                    cont_sg_norm >= 0 ~ "p"),
           cont_sg_sig = case_when((cont_sg_norm_lci < 0 & 0 < cont_sg_norm_uci) ~ "N",
                                   TRUE ~ "Y"),
           cont_sg_display = case_when(cont_sg_sig == "Y" ~ cont_sg_norm,
                                       T ~ NA_real_),
           cont_sg_code = factor(case_when(
             cont_sg_sig == "Y" & cont_sg_sign == "n" ~ "Neg",
             cont_sg_sig == "Y" & cont_sg_sign == "p" ~ "Pos",
             cont_sg_sig == "N" ~ "Non Sig"), 
             levels=c("Neg", "Pos", "Non Sig")),
           
           cont_ghm_sign = case_when(cont_ghm_scale < 0 ~ "n",
                                     cont_ghm_scale >= 0 ~ "p"),
           cont_ghm_sig = case_when((cont_ghm_scale_lci < 0 & 0 < cont_ghm_scale_uci) ~ "N",
                                    TRUE ~ "Y"),
           cont_ghm_display = case_when(cont_ghm_sig == "Y" ~ cont_ghm_scale,
                                        T ~ NA_real_),
           cont_ghm_code = factor(case_when(cont_ghm_sig == "Y" & cont_ghm_sign == "n" ~ "Neg",
                                            cont_ghm_sig == "Y" & cont_ghm_sign == "p" ~ "Pos",
                                            cont_ghm_sig == "N" ~ "Non Sig"), 
                                  levels=c("Neg", "Pos", "Non Sig")))
  
  
  # INTERACTION MODEL
  load(int_modlist_full[i])
  fe <- fixef(out$model)
  intdf <- tibble("species" = out$species,
                  "inter" = fe["sg_norm:ghm_scale", "Estimate"],
                  "inter_lci" = fe["sg_norm:ghm_scale", "Q2.5"],
                  "inter_uci" = fe["sg_norm:ghm_scale", "Q97.5"]) %>%
    mutate(inter_sign = case_when(inter < 0 ~ "n",
                                  inter >= 0 ~ "p"),
           inter_sig = case_when((inter_lci < 0 & 0 < inter_uci) ~ "N",
                                 TRUE ~ "Y"),
           inter_display = case_when(inter_sig == "Y" ~ inter,
                                     T ~ NA_real_),
           int_code = factor(case_when(inter_sig == "Y" & inter_sign == "n" ~ "Neg",
                                       inter_sig == "Y" & inter_sign == "p" ~ "Pos",
                                       inter_sig == "N" ~ "Non Sig"),
                             levels=c("Neg", "Pos", "Non Sig")))
  # get observed quantiles of ghm to set "low" and "high" human mod
  ghmq <- quantile(out$data$ghm_scale, probs = c(0.10, 0.90), na.rm = T)
  
  # Conditional Effects Plot for interaction
  ce_int <- conditional_effects(x=out$model,
                                effects = "sg_norm:ghm_scale",
                                int_conditions = list(ghm_scale = ghmq),
                                re_formula = NA)
  
  if(intdf$inter_sig == "Y"){ # If the effect is significant...
    (int_ce_plot <-  plot(ce_int, plot = FALSE,
                          line_args = list("se"=F,
                                           "size" = 8))[[1]] +
       theme_tufte() +
       scale_color_manual(values = palnew, name = "Human \n Modification",
                          labels = c("High", "Low")) +
       scale_fill_manual(values = palnew, name = "Human \n Modification",
                         labels = c("High", "Low")) +
       theme(axis.line = element_line(size = 8),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             aspect.ratio = 1,
             legend.position = "none"
             
       ))
    ggsave(filename = glue("out/niche_plots/{out$species}_niche_int.png"), int_ce_plot)
    plint <- 1
  }else{ # ...if the effect is not significant
    (int_ce_plot <-  plot(ce_int, plot = FALSE,
                          line_args = list("se"=F, "size" = 8))[[1]] +
       theme_tufte() +
       # xlab("Human Mobility") +
       # ylab("")+
       scale_color_grey(name = "Human \n Modification",
                        labels = c("High", "Low"))+
       scale_fill_grey(name = "Human \n Modification",
                       labels = c("High", "Low"))+
       theme(axis.line = element_line(size = 8),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             aspect.ratio = 1,
             legend.position = "none"))
    
    plint <- 0
  }
  
  # gather results tables
  res_out[[i]] <- adddf %>% 
    left_join(intdf) %>% 
    left_join(contdf)
  
  
} # i

# combine dfs
res_out_df <- do.call("rbind", res_out)
write_csv(x=res_out_df, file = glue("out/niche_mod_summary_{Sys.Date()}-noenv.csv"))

# combine plots and save out
row_reduced <- row[which(pl==1)]
fin_plots <- do.call(c, row_reduced)

(comb <- wrap_plots(fin_plots, byrow = T, ncol = 4, guides = "collect", 
                    widths=c(3,5,5,5), heights = c(5)))
ggsave(filename='out/niche_plots/comb_08092022.pdf', width = 7, height = 11.5, 
       plot = comb)



#---- Trait Model ----#

# Load trait model
load("out/single_species_models/area_trait/size_trait_mod_2022-07-30.rdata")

conds <- make_conditions(out$model, vars = "diet")

ghmq_trait <- quantile(out$data$ghm_scale, probs = c(0.05, 0.95), na.rm = T)

sgq_trait <- quantile(out$data$sg_norm, probs = c(0.10, 0.90), na.rm = T)

ce_trait <- conditional_effects(x=out$model, 
                                effects = "sg_scale:ghm_scale",
                                int_conditions = list(ghm_scale = ghmq_trait
                                                      # sg_norm = c(-.2225,-0.2138)
                                ),
                                conditions = conds,
                                re_formula = NA)

plot(ce_trait,
     facet_args = list("scales" = "free_y"))[[1]]+
  # xlim(c(-.2225,1))+
  theme_tufte()+
  theme(axis.line = element_line(size = .5),
        # axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

ce_trait %>% 
  mutate(diet_sig = case_when(diet = ))

#---- Save output ---#
# message(glue('Saving to {.outPF}'))

message(glue('Script complete in {diffmin(t0)} minutes'))


#---- SCRATCH ----#

