'
Plot Effect Sizes.  Generate figure showing effects sizes for both 

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
  
  .nichePF <- file.path(.wd,'out/niche_mod_summary_2022-08-26.csv')
  .areaPF <- file.path(.wd,'out/area_mod_summary_2022-08-26.csv')
  .traitPF <- file.path(.wd, 'raw_data/anthropause_data_sheet.csv')
  .outPF <- file.path(.wd,'out/')
  
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

pal <- c("#D3D3D3", "#5D3A9B")


#---- Load Data ----#

traits <- read_csv(.traitPF)
niche0 <- read_csv(.nichePF)
# area0 <- read_csv(.areaPF)


#---- Process Data ----#

niche <- niche0 %>% 
  filter(species != "Anser") %>% 
  left_join(traits, by = c("species" = "Species")) %>% 
  mutate(species_ghm = fct_reorder(species, ghm_scale, min),
         species_sg = fct_reorder(species, sg_norm, min)) 

# Make color vector for x axis labels
ghm_sort <- niche %>% 
  arrange(ghm_scale) # sort df by the ghm vals
niche_ghm_labs <- ifelse(ghm_sort$ghm_sig == "N", pal[1], pal[2]) # make color vector

# Make color vector for x axis labels
sg_sort <- niche %>% 
  arrange(sg_norm) # sort df by the ghm vals
niche_sg_labs <- ifelse(sg_sort$sg_sig == "N", pal[1], pal[2]) # make color vector

#---- Make Plots ----#

#-- Niches --#

# GHM
(niche_ghm_plot <- ggplot(niche)+
  geom_pointrange(aes(y = ghm_scale, ymin = ghm_scale_lci, ymax = ghm_scale_uci, 
                      x= species_ghm, color = ghm_sig)) +
  scale_color_manual(values = pal)+
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  xlab("Species")+
  ylab("Standardized Effect Size") +
  ggtitle("Effects of Human Modification on Niche Breadth")+
  theme_tufte()+
  theme(axis.text.x=element_text(angle=75, hjust = 1, colour = niche_ghm_labs),
        legend.position = "none")+
# facet_wrap(~species_group)
NULL)

# SG
(niche_sg_plot <- ggplot(niche)+
    geom_pointrange(aes(y = sg_norm, ymin = sg_norm_lci, ymax = sg_norm_uci, 
                        x= species_sg, color = sg_sig)) +
    scale_color_manual(values = pal)+
    geom_hline(aes(yintercept = 0), linetype = "dashed")+
    xlab("Species")+
    ylab("Standardized Effect Size") +
    ggtitle("Effects of Human Mobility on Niche Breadth")+
    theme_tufte()+
    theme(axis.text.x=element_text(angle=75, hjust = 1, colour = niche_sg_labs),
          legend.position = "none")+
    # facet_wrap(~species_group)
    NULL)

niche_ghm_plot+niche_sg_plot
