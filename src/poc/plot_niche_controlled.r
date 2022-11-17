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


message('Loading  models...')
cont_modlist <- list.files(path=file.path(.datP, "niche_controlled/"), full.names = T)

#---- Make Plots ----#

pal2 <- c("#E66100", "#5D3A9B") # 2 color pallete
pal3 <- c(pal2, "#808080") # add gray to the pallete
palnew <- c("#7552A3", "#CEBEDA")
palgray <- c("#808080", "#D3D3D3")




#---- Gather Results ----#

res_out <- list()
row <- list()
pl <- c()
for(i in 1:length(modlist)){
  tmp <- list()
  
  # CONTROLLED MODEL
  load(cont_modlist[i]) # load model
  fe <- fixef(out$model) #get fixed effects
  re <- posterior_summary(out$model, variable = c("sd_grp__Intercept", "sigma"))
  df <- tibble("species"=out$species, # grab estimates
               
               # AREA EFFECTS
               "area"=as.numeric(fe["log_area_scale", "Estimate"]),
               "area_lci"=fe["log_area_scale", "Q2.5"],
               "area_uci"=fe["log_area_scale", "Q97.5"],
               
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
    mutate(area_sign = case_when(area < 0 ~ "n",
                               area >= 0 ~ "p"),
           area_sig = case_when((area_lci < 0 & 0 < area_uci) ~ "N",
                              TRUE ~ "Y"),
           area_display = case_when(area_sig == "Y" ~ area,
                                  T ~ NA_real_),
           code = factor(case_when(area_sig == "Y" & area_sign == "n" ~ "Neg",
                                   area_sig == "Y" & area_sign == "p" ~ "Pos",
                                   area_sig == "N" ~ "Non Sig"), 
                         levels=c("Neg", "Pos", "Non Sig")),
           
           sg_sign = case_when(sg_norm < 0 ~ "n",
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
                         levels=c("Neg", "Pos", "Non Sig")),
           
           ICC = group/(group + resid),
           var_ratio = resid/group)
  
  # gather results tables
  res_out[[i]] <- df
  
}

# combine dfs
res_out_df <- do.call("rbind", res_out)
write_csv(x=res_out_df, file = glue("out/niche_controlled_summary_{Sys.Date()}.csv"))



#---- Plots ----#

# Define color palette
pal <- c("#D3D3D3", "#5D3A9B")

# Make vector of species names for plotting
sp_vec <- c("Puma concolor", "Canis lupus", "Canis latrans", "Lynx rufus",
            "Ursus americanus", 
            
            "Alces alces", "Antilocapra americana",
            "Cervus elaphus", 
            "Odocoileus hemionus",
            "Odocoileus virginianus",
            
            "Anas acuta", "Anas clypeata", "Anas cyanoptera", "Anas strepera",
            "Anas americana", "Anas platyrhynchos", "Anas crecca",
            "Anser albifrons", "Chen rossii", "Chen caerulescens", 
            "Grus canadensis", 
            
            "Circus cyaneus", "Corvus corax", "Haliaeetus leucocephalus")

dat <- res_out_df %>% 
  filter(species != "Anser") %>% 
  # left_join(traits, by = c("species" = "Species")) %>% 
  mutate(species_ghm = fct_reorder(species, ghm_scale, min),
         species_sg = fct_reorder(species, sg_norm, min),
         species_man = fct_relevel(species, sp_vec),
         group = case_when(species == "Puma concolor" ~ "A", 
                           species == "Canis lupus" ~ "A", 
                           species == "Canis latrans" ~ "A", 
                           species == "Lynx rufus" ~ "A",
                           species == "Ursus americanus" ~ "A", 
                           
                           species == "Alces alces" ~ "B", 
                           species == "Antilocapra americana" ~ "B",
                           species == "Cervus elaphus" ~ "B", 
                           species == "Odocoileus hemionus" ~ "B",
                           species == "Odocoileus virginianus" ~ "B",
                           
                           species == "Anas acuta" ~ "C", 
                           species == "Anas clypeata" ~ "C", 
                           species == "Anas cyanoptera" ~ "C", 
                           species == "Anas strepera" ~ "C",
                           species == "Anas americana" ~ "C", 
                           species == "Anas platyrhynchos" ~ "C", 
                           species == "Anas crecca" ~ "C",
                           species == "Anser albifrons" ~ "C", 
                           species == "Chen rossii" ~ "C", 
                           species == "Chen caerulescens" ~ "C", 
                           species == "Grus canadensis" ~ "C", 
                           
                           species == "Circus cyaneus" ~ "D", 
                           species == "Corvus corax" ~ "D", 
                           species == "Haliaeetus leucocephalus" ~ "D")) 

# Space Use
(cont_area_plot <- ggplot(dat)+
    geom_pointrange(aes(x = area, xmin = area_lci, xmax = area_uci, 
                        y= species_man, group = species_man, color = area_sig)) +
    scale_color_manual(values = pal)+
    geom_vline(aes(xintercept = 0), linetype = "dashed")+
    ylab("Species")+
    xlab("Standardized Effect Size") +
    ggtitle("Space Use")+
    scale_y_discrete(limits=rev) +
    # make space to accommodate the fake axis
    expand_limits(x = -3.1) +
    # create a strip of white background under the fake axis
    geom_rect(xmin = -1.8, xmax = -2.8, ymin = 0, ymax = nrow(dat) + 1, fill = "white") +
    # fake axis layer, aligned below y = 0
    geom_text(aes(label = species, y = species_man, color = area_sig), x = -0.5, hjust = 1) +
    # specify the font colours for fake axis
    # scale_colour_manual(values = c("black", "red"), guide = F) +
    # hide the actual x-axis text / tick
    facet_grid(rows = vars(group), space = "free_y", scales = "free_y") +    
    theme_tufte()+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          # axis.text.y=element_text(colour = rev(niche_ghm_labs_man)),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text.y = element_blank())+
    NULL)

# GHM
(cont_ghm_plot <- ggplot(dat)+
    geom_pointrange(aes(x = ghm_scale, xmin = ghm_scale_lci, xmax = ghm_scale_uci, 
                        y= species_man, group = species_man, color = ghm_sig)) +
    scale_color_manual(values = pal)+
    geom_vline(aes(xintercept = 0), linetype = "dashed")+
    ylab("Species")+
    xlab("Standardized Effect Size") +
    ggtitle("Modification")+
    scale_y_discrete(limits=rev) +
    # make space to accommodate the fake axis
    expand_limits(x = -2.75) +
    # create a strip of white background under the fake axis
    geom_rect(xmin = -1, xmax = -2, ymin = 0, ymax = nrow(dat) + 1, fill = "white") +
    # fake axis layer, aligned below y = 0
    geom_text(aes(label = species, y = species_man, color = ghm_sig), x = -1, hjust = 1) +
    # specify the font colours for fake axis
    # scale_colour_manual(values = c("black", "red"), guide = F) +
    # hide the actual x-axis text / tick
    facet_grid(rows = vars(group), space = "free_y", scales = "free_y") +    
    theme_tufte()+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          # axis.text.y=element_text(colour = rev(niche_ghm_labs_man)),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text.y = element_blank())+
    NULL)



# SG
(cont_sg_plot <- ggplot(dat)+
    geom_pointrange(aes(x = sg_norm, xmin = sg_norm_lci, xmax = sg_norm_uci, 
                        y= species_man, group = species_man, color = sg_sig)) +
    scale_color_manual(values = pal)+
    geom_vline(aes(xintercept = 0), linetype = "dashed")+
    ylab("Species")+
    xlab("Standardized Effect Size") +
    ggtitle("Mobility")+
    scale_y_discrete(limits=rev) +
    # make space to accommodate the fake axis
    expand_limits(x = -2.75) +
    # create a strip of white background under the fake axis
    geom_rect(xmin = -1.8, xmax = -2.5, ymin = 0, ymax = nrow(dat) + 1, fill = "white") +
    # fake axis layer, aligned below y = 0
    geom_text(aes(label = species, y = species_man, color = sg_sig), x = -1.8, hjust = 1) +
    # specify the font colours for fake axis
    # scale_colour_manual(values = c("black", "red"), guide = F) +
    # hide the actual x-axis text / tick
    facet_grid(rows = vars(group), space = "free_y", scales = "free_y") +    
    theme_tufte()+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          # axis.text.y=element_text(colour = rev(niche_ghm_labs_man)),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text.y = element_blank())+
    NULL)


#-- Combine and Save --#

(niche_control_plot <- cont_area_plot + cont_ghm_plot + cont_sg_plot + plot_annotation(
  title = 'Niche Breadth'))

ggsave(plot = niche_control_plot, filename = glue("figs/niche_control_plot_{Sys.Date()}.png"), 
       width = 20, height = 5)



#-- SCRATCH --#

ggplot(dat) +
  # geom_pointrange(aes(x = ICC, y = area, ymin = area_lci, ymax = area_uci))
  geom_point(aes(x = ICC, y = ghm_scale))
