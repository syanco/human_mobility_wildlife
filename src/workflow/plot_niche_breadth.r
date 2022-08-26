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

# Interaction Models
message('Loading interaction models...')
int_modlist <- list.files( path=file.path(.datP, "niche/"), full.names=TRUE ) 


# SG Models
message('Loading safegraph models...')

sg_modlist <- list.files( path=file.path(.datP, "niche_sg/"), full.names=TRUE ) 


# SG Models
message('Loading safegraph models...')

ghm_modlist <- list.files( path=file.path(.datP, "niche_ghm/"), full.names=TRUE ) 

#---- Make Plots ----#

# Define some color palettes
pal2 <- c("#E66100", "#5D3A9B") # 2 color palette
pal3 <- c(pal2, "#808080") # add gray to the palette
palnew <- palnew <- c("#7552A3", "#CEBEDA")
palgray <- c("#808080", "#D3D3D3")


#-- make combined plots --#

# !!!!!!!!! #
# TODO:  this just manually matches the species lists... 
# Need to make a better way to ensure I'm calling the same species across all
# model types - maybe use species list from DB and then call model files by 
# regex?
# !!!!!!!!! #

# sg_modlist <- sg_modlist[-22]

# !!!!!!!!! !

res_out <- list()
row <- list()
pl <- c()
for(i in 1:length(sg_modlist)){
  tmp <- list()
  
  # SG ONLY
  load(sg_modlist[i]) # load model
  fe <- fixef(out$model) #get fixed effects
  sgdf <- tibble("species"=out$species, # grab estimates
                 "sg_norm"=as.numeric(fe["sg_norm", "Estimate"]),
                 "sg_norm_lci"=fe["sg_norm", "Q2.5"],
                 "sg_norm_uci"=fe["sg_norm", "Q97.5"]) %>% 
    mutate(sg_sign = case_when(sg_norm < 0 ~ "n",
                               sg_norm >= 0 ~ "p"),
           sg_sig = case_when((sg_norm_lci < 0 & 0 < sg_norm_uci) ~ "N",
                              TRUE ~ "Y"),
           sg_display = case_when(sg_sig == "Y" ~ sg_norm,
                                  T ~ NA_real_),
           sg_code = factor(case_when(sg_sig == "Y" & sg_sign == "n" ~ "Neg",
                                   sg_sig == "Y" & sg_sign == "p" ~ "Pos",
                                   sg_sig == "N" ~ "Non Sig"), 
                         levels=c("Neg", "Pos", "Non Sig")))
  
  # Get conditional effects
  ce_sg <- conditional_effects(x=out$model, 
                               effects = "sg_norm",
                               re_formula = NA)
  # if(sgdf$sg_sig == "N"){next}
  if(sgdf$sg_sig == "Y"){ # If the effect is significant...
    (sg_ce_plot <-  plot(ce_sg, plot = F,
                         line_args = list("se" = F,
                                          "color" = pal2[2]))[[1]] + 
       # scale_color_manual(values = palnew[3])+         
       theme_tufte() +
       xlab("Human Mobility") +
       ylab("Space Use")+
       theme(axis.line = element_line(size = .5),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             # aspect.ratio = 1
             ))
    plsg <- 1
  }else{ # ...not significant:
    (sg_ce_plot <-  plot(ce_sg, plot = F,
                         line_args = list("se" = F,
                                          "color" = "#808080"))[[1]] + 
       # scale_color_manual(values = palnew[3])+         
       theme_tufte() +
       xlab("Human Mobility") +
       ylab("Space Use")+
       theme(axis.line = element_line(size = .5),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             # aspect.ratio = 1
             ))
    plsg <- 0
  }
  
  
  # GHM ONLY
  load(ghm_modlist[i]) # load model
  fe <- fixef(out$model) #get fixed effects
  ghmdf <- tibble("species"=out$species, # grab estimates
                  "ghm_scale"=as.numeric(fe["ghm_scale", "Estimate"]),
                  "ghm_scale_lci"=fe["ghm_scale", "Q2.5"],
                  "ghm_scale_uci"=fe["ghm_scale", "Q97.5"]) %>% 
    mutate(ghm_sign = case_when(ghm_scale < 0 ~ "n",
                                ghm_scale >= 0 ~ "p"),
           ghm_sig = case_when((ghm_scale_lci < 0 & 0 < ghm_scale_uci) ~ "N",
                               TRUE ~ "Y"),
           ghm_display = case_when(ghm_sig == "Y" ~ ghm_scale,
                                   T ~ NA_real_),
           ghm_code = factor(case_when(ghm_sig == "Y" & ghm_sign == "n" ~ "Neg",
                                   ghm_sig == "Y" & ghm_sign == "p" ~ "Pos",
                                   ghm_sig == "N" ~ "Non Sig"), 
                         levels=c("Neg", "Pos", "Non Sig")))
  
  # Get conditional effects
  ce_ghm <- conditional_effects(x=out$model, 
                                effects = "ghm_scale",
                                re_formula = NA)
  # if(sgdf$sg_sig == "N"){next}
  if(ghmdf$ghm_sig == "Y"){ # If the effect is significant...
    (ghm_ce_plot <-  plot(ce_ghm, plot = F,
                          line_args = list("se" = F,
                                           "color" = pal2[2]))[[1]] + 
       # scale_color_manual(values = palnew[3])+         
       theme_tufte() +
       xlab("Human Modification") +
       ylab("Space Use")+
       theme(axis.line = element_line(size = .5),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             # aspect.ratio = 1
       ))
    plghm <- 1
  }else{ # ...not significant:
    (ghm_ce_plot <-  plot(ce_ghm, plot = F,
                          line_args = list("se" = F,
                                           "color" = "#808080"))[[1]] + 
       # scale_color_manual(values = palnew[3])+         
       theme_tufte() +
       xlab("Human Modification") +
       ylab("Space Use")+
       theme(axis.line = element_line(size = .5),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             # aspect.ratio = 1
       ))
    plghm <- 0
  }
  
  # INTERACTION MODEL
  load(int_modlist[i])
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
                                           "size" = 4))[[1]] +
       theme_tufte() +
       scale_color_manual(values = palnew, name = "Human \n Modification",
                          labels = c("High", "Low")) +
       scale_fill_manual(values = palnew, name = "Human \n Modification",
                         labels = c("High", "Low")) +
       theme(axis.line = element_line(size = 4),
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
                          line_args = list("se"=F, "size" = 4))[[1]] +
       theme_tufte() +
       # xlab("Human Mobility") +
       # ylab("")+
       scale_color_grey(name = "Human \n Modification",
                        labels = c("High", "Low"))+
       scale_fill_grey(name = "Human \n Modification",
                       labels = c("High", "Low"))+
       theme(axis.line = element_line(size = 4),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             axis.title = element_blank(),
             aspect.ratio = 1,
             legend.position = "none"))

    plint <- 0
  }
  
  # gather results tables
  res_out[[i]] <- ghmdf %>% left_join(sgdf) %>% left_join(intdf)
  
  # gather plot objects
  tmp[[1]] <- wrap_elements(textGrob(glue("{out$species}"),
                                     gp = gpar(fontsize = 6)))
  tmp[[2]] <- sg_ce_plot 
  tmp[[3]] <- ghm_ce_plot
  tmp[[4]] <- int_ce_plot
  
  # +
  # plot_layout(widths=c(1,1,1), heights = c(1))
  row[[i]] <- tmp 
  # wrap_plots(tmp, widths=c(1,1,1), heights = c(1))
  if(plsg == 1 | plghm == 1 | plint ==1){pl[i] <- 1}else(pl[i] <- 0)
  
} # i

# combine dfs
res_out_df <- do.call("rbind", res_out)
write_csv(x=res_out_df, file = glue(""))

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

