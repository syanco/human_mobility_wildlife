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
  .areaPF <- file.path(.wd,'out/area_mod_summary_2022-08-31.csv')
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
area0 <- read_csv(.areaPF)


#---- Process Data ----#

# Make vector of species names for plotting

sp_vec <- c("Puma concolor", "Canis lupus", "Canis latrans", "Lynx rufus",
            "Ursus americanus", 
            
            "Alces alces", "Antilocapra americana",
            "Cervus elaphus", "Odocoileus hemionus", "Odocoileus virginianus",
            
            "Anas acuta", "Anas clypeata", "Anas cyanoptera", "Anas strepera",
            "Anas americana", "Anas platyrhynchos", "Anas crecca",
            "Anser albifrons", "Chen rossii", "Chen caerulescens", 
            "Grus canadensis", 
            
            "Circus cyaneus", "Corvus corax", "Haliaeetus leucocephalus")

niche <- niche0 %>% 
  filter(species != "Anser") %>% 
  left_join(traits, by = c("species" = "Species")) %>% 
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

area <- area0 %>% 
  # filter(species != "Anser") %>% 
  left_join(traits, by = c("species" = "Species")) %>% 
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




#---- Make Plots ----#

#-- Niches --#
# Make color vectors for x axis labels
# by GHM est
ghm_sort <- niche %>% 
  arrange(ghm_scale) # sort df by the ghm vals
niche_ghm_labs <- ifelse(ghm_sort$ghm_sig == "N", pal[1], pal[2]) # make color vector

# by sg est
sg_sort <- niche %>% 
  arrange(sg_norm) # sort df by the ghm vals
niche_sg_labs <- ifelse(sg_sort$sg_sig == "N", pal[1], pal[2]) # make color vector

# manual sort *match fig 1
man_sort <- niche %>% 
  arrange(group, species_man) # sort df by the factor levels
niche_sg_labs_man <- ifelse(man_sort$sg_sig == "N", pal[1], pal[2]) # make color vector
niche_ghm_labs_man <- ifelse(man_sort$ghm_sig == "N", pal[1], pal[2]) # make color vector

# GHM
(niche_ghm_plot <- ggplot(niche)+
    geom_pointrange(aes(x = ghm_scale, xmin = ghm_scale_lci, xmax = ghm_scale_uci, 
                        y= species_man, group = species_man, color = ghm_sig)) +
    scale_color_manual(values = pal)+
    geom_vline(aes(xintercept = 0), linetype = "dashed")+
    ylab("Species")+
    xlab("Standardized Effect Size") +
    ggtitle("Human Modification")+
    scale_y_discrete(limits=rev) +
    # make space to accommodate the fake axis
    expand_limits(x = -1.6) +
    # create a strip of white background under the fake axis
    geom_rect(xmin = -1.5, xmax = -2, ymin = 0, ymax = nrow(niche) + 1, fill = "white") +
    # fake axis layer, aligned below y = 0
    geom_text(aes(label = species, y = species_man, color = ghm_sig), x = -1.4, hjust = 1) +
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
(niche_sg_plot <- ggplot(niche)+
    geom_pointrange(aes(x = sg_norm, xmin = sg_norm_lci, xmax = sg_norm_uci, 
                        y= species_man, group = species_man, color = sg_sig)) +
    scale_color_manual(values = pal)+
    geom_vline(aes(xintercept = 0), linetype = "dashed")+
    ylab("Species")+
    xlab("Standardized Effect Size") +
    ggtitle("Human Mobility")+
    scale_y_discrete(limits=rev) +
    # make space to accommodate the fake axis
    expand_limits(x = -2.6) +
    # create a strip of white background under the fake axis
    geom_rect(xmin = -1.5, xmax = -2, ymin = 0, ymax = nrow(niche) + 1, fill = "white") +
    # fake axis layer, aligned below y = 0
    geom_text(aes(label = species, y = species_man, color = sg_sig), x = -3.5, hjust = 1) +
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

niche_effect_plot <- niche_ghm_plot+niche_sg_plot + plot_annotation(
  title = 'Anthropogenic impacts on animal niche breadth')
ggsave(plot = niche_effect_plot, filename = glue("figs/niche_effect_size_plot_{Sys.Date()}.png"), 
       width = 20, height = 5)


#-- Area --#

# Make color vectors for x axis labels
# by GHM est
ghm_sort_area <- area %>% 
  arrange(ghm_scale) # sort df by the ghm vals
area_ghm_labs <- ifelse(ghm_sort_area$ghm_sig == "N", pal[1], pal[2]) # make color vector

# by sg est
sg_sort_area <- area %>% 
  arrange(sg_norm) # sort df by the ghm vals
area_sg_labs <- ifelse(sg_sort_area$sg_sig == "N", pal[1], pal[2]) # make color vector

# manual sort *match fig 1
man_sort_area <- area %>% 
  arrange(group, species_man) # sort df by the factor levels
area_sg_labs_man <- ifelse(man_sort_area$sg_sig == "N", pal[1], pal[2]) # make color vector
area_ghm_labs_man <- ifelse(man_sort_area$ghm_sig == "N", pal[1], pal[2]) # make color vector

# GHM
(area_ghm_plot <- ggplot(area)+
    geom_pointrange(aes(x = ghm_scale, xmin = ghm_scale_lci, xmax = ghm_scale_uci, 
                        y= species_man, group = species_man, color = ghm_sig)) +
    scale_color_manual(values = pal)+
    geom_vline(aes(xintercept = 0), linetype = "dashed")+
    ylab("Species")+
    xlab("Standardized Effect Size") +
    ggtitle("Human Modification")+
    scale_y_discrete(limits=rev) +
    # make space to accommodate the fake axis
    expand_limits(x = -1.4) +
    # create a strip of white background under the fake axis
    geom_rect(xmin = -1.5, xmax = -2, ymin = 0, ymax = nrow(niche) + 1, fill = "white") +
    # fake axis layer, aligned below y = 0
    geom_text(aes(label = species, y = species_man, color = ghm_sig), x = -0.6, hjust = 1) +
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
(area_sg_plot <- ggplot(area)+
    geom_pointrange(aes(x = sg_norm, xmin = sg_norm_lci, xmax = sg_norm_uci, 
                        y= species_man, group = species_man, color = sg_sig)) +
    scale_color_manual(values = pal)+
    geom_vline(aes(xintercept = 0), linetype = "dashed")+
    ylab("Species")+
    xlab("Standardized Effect Size") +
    ggtitle("Human Mobility")+
    scale_y_discrete(limits=rev) +
    # make space to accommodate the fake axis
    expand_limits(x = -1) +
    # create a strip of white background under the fake axis
    geom_rect(xmin = -1.5, xmax = -2, ymin = 0, ymax = nrow(niche) + 1, fill = "white") +
    # fake axis layer, aligned below y = 0
    geom_text(aes(label = species, y = species_man, color = sg_sig), x = -0.5, hjust = 1) +
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

area_effect_plot <- area_ghm_plot+area_sg_plot+ plot_annotation(
  title = 'Anthropogenic impacts on animal space use')
ggsave(plot = area_effect_plot, filename = glue("figs/area_effect_size_plot_{Sys.Date()}.png"), 
       width = 20, height = 5)

