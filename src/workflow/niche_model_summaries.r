#!/usr/bin/env Rscript 
# Plot model summary sheet


#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- getwd()
  
  .datP <- file.path(.wd,'out/single_species_models')
  .outPF <- file.path(.wd,'figs/area_fig.png')
  .dbPF <- file.path(.wd,'processed_data/mosey_mod_2023.db')
  
} else {
  library(docopt)
  # library(rprojroot)
  
  # ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  
  source('src/funs/input_parse.r')
  
  #.list <- trimws(unlist(strsplit(ag$list,',')))
  # .datP <- makePath(ag$dat)
  .datP <- file.path(.wd,'out/single_species_models')
  # .outPF <- makePath(ag$out)
  .dbPF <- '/tmp/mosey_mod.db'
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
    library(gridExtra)
    library(rnaturalearth)
    library(DBI)
    library(RSQLite)
    library(sf)
    library(ggsflabel)
    library(janitor)
    library(glue)
  }))

#Source all files in the auto load funs directory
list.files('src/funs/auto',full.names=TRUE) %>%
  walk(source)

palnew <- c("#F98177", "#8895BF")

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

# Get US Background
us <- ne_states(country = "United States of America", returnclass = "sf") #%>%
        #st_set_crs(4326)


#---- Initialize database ----#
message("Connecting to database...")
invisible(assert_that(file.exists(.dbPF)))

db <- dbConnect(RSQLite::SQLite(), .dbPF)

invisible(assert_that(length(dbListTables(db))>0))

#---- Perform analysis ----#

# dbListTables(db)

#-- Load Cleaned Data
evt0 <- tbl(db, 'event_clean') %>% collect()
std0 <- tbl(db, 'study_clean') %>% 
  collect() 
# beepr::beep()

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
        
        #- Model Basics -#        
        load(add_modlist_full[i]) # load model
        addmod <- out$model
        sp <- out$species
        out_add <- out
        fe_add <- fixef(out_add$model) #get fixed effects
        mod <- "Additive Model"
        
        data_sum_tbl <- out_add$data %>% 
          group_by(studyid) %>% 
          summarize(num_inds = n_distinct(ind_f),
                    #num_weeks = n_distinct(wk),
                    num_obs = n()) %>%
          mutate(obs_per_ind = num_obs/num_inds) %>% 
          adorn_totals() %>% 
          tableGrob()
        
        
        #- Make conditional predictions -#
        
        # get observed quantiles of ghm to set "low" and "high" human mod
        ghmq <- quantile(out$data$ghm_scale, probs = c(0.1, 0.9), na.rm = T)
        sgq <- quantile(out$data$sg_norm, probs = c(0.1, 0.9), na.rm = T)
        
        # Conditional Effects Plot for safegraph
        ce_add_sg <- conditional_effects(x=addmod, 
                                         effects = c("sg_norm"),
                                         int_conditions = list(ghm_scale = ghmq,
                                                               sg_norm = sgq),
                                         re_formula = NA, 
                                         plot = F) 
        # Plot
        (add_ce_plot_sg <-  plot(ce_add_sg, 
                                 plot = F,
                                 rug = F,
                                 line_args = list("se" = T))[[1]] + 
            xlab("Human Mobility") +
            ylab("Scaled Niche Breadth")+
            theme_minimal() +
            theme(panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  # legend.position = "none",
                  axis.title = element_text(size = 11,
                                            face = "bold"),
                  axis.ticks = element_line(color = "#4a4e4d"),
                  axis.text = element_text(size = 12),
                  text = element_text(family = "Arial", color = "#4a4e4d")) 
        )
        
        # Conditional Effects Plot for GHM
        ce_add_ghm <- conditional_effects(x=addmod, 
                                          effects = c("ghm_scale"),
                                          int_conditions = list(ghm_scale = ghmq,
                                                                sg_norm = sgq),
                                          re_formula = NA, 
                                          plot = F) 
        # Plot
        (add_ce_plot_ghm <-  plot(ce_add_ghm, 
                                  plot = F,
                                  rug = F,
                                  line_args = list("se" = T))[[1]] + 
            xlab("Human Modification") +
            ylab("Scaled Niche Breadth")+
            theme_minimal() +
            theme(panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  # legend.position = "none",
                  axis.title = element_text(size = 11,
                                            face = "bold"),
                  axis.ticks = element_line(color = "#4a4e4d"),
                  axis.text = element_text(size = 12),
                  text = element_text(family = "Arial", color = "#4a4e4d")) 
        )
        
        #- Extract parameter table as grob -#
        param_table <- tableGrob(parameters(addmod))
        
        #- Posterior Predictive Plot -#
        (pp_dens <- pp_check(addmod)+ggtitle("Posterior Predictive Distribution"))
        
        #- Predictive Error Plot -#
        (pp_err <- pp_check(addmod, type='error_scatter_avg')+ggtitle("Posterior Predictive Errors"))  
        
        #- MCMC Trace Plot -#
        (trace <- mcmc_plot(addmod, type = "trace") +ggtitle("MCMC Traces"))    
        
        
        
        #-- Maps --#
        
        #- Daily Fixes Map -#
        
        # Filter to sp
        evt_sp <- evt0  %>% 
          filter(species==sp) %>% 
          mutate(doy = day(timestamp)) 
        
        # Summarize to daily fix, covert to sf
        evt_daily <- evt_sp %>% 
          st_as_sf(coords = c("lon", "lat"), crs=4326) %>% 
          group_by(study_id, ind_f, doy) %>% 
          summarize(st_union(geometry),
                    study_id = study_id[1]) %>% 
          st_centroid()
        
        # Define bounding range
        xrange <- c(min(evt_sp$lon)-1, max(evt_sp$lon)+1)
        yrange <- c(min(evt_sp$lat)-1, max(evt_sp$lat)+1)
        
        # Plot map
        (zoom_map <- ggplot()+
            geom_sf(data = us)+
            geom_sf(data = evt_daily, inherit.aes = F, aes(color = ind_f))+
            # geom_sf_label_repel(data = std_sf, aes(label = species), force_pull = 0, force = 10) +
            coord_sf() + 
            xlim(xrange) +
            ylim(yrange)+
            ggtitle("Individual Daily Positions")+
            theme(legend.position = "none") )
        
        #- Study Map -#
        
        # Get study centroids (based on points)
        std_sf <- evt_sp %>% 
          left_join(std0) %>% 
          st_as_sf(coords = c("lon", "lat"), crs=4326) %>% 
          group_by(study_id) %>% 
          summarize(st_union(geometry),
                    study_id = study_id[1],
                    study_name = study_name[1]) %>% 
          st_centroid()
        
        # Get bounding box of centroids to set plot lims
        stdbb <- st_bbox(std_sf)

        # ensure CRS match before plotting
        CRS_same <- st_crs(us) == st_crs(std_sf)
        if (CRS_same) {
          message("CRS of std_sf and us is the same.")
        } else {
          us <- st_transform(us, st_crs(std_sf))
          message("Transformed CRS of us polygons to CRS of std_sf")
        }
        
        # make plot
        study_plot <- ggplot()+
          geom_sf(data = us)+
          geom_sf(data = std_sf, inherit.aes = F)+
          xlim(c(stdbb['xmin']-1, stdbb['xmax']+1)) +
          ylim(c(stdbb['ymin']-1, stdbb['ymax']+1))+
          ggtitle("Study Centroids") +
          geom_sf_label_repel(data = std_sf, aes(label = study_name), force_pull = 0, force = 10) +
          coord_sf()
        
        
        #---- Assemble Plots ---#
        
        # Design layout
        design <- " #IIII#
                    #IIII#
                    AAABBB
                    AAABBB
                    AAABBB
                    CCCCCC
                    CCCCCC
                    EEEFFF
                    EEEFFF
                    GGGGGG
                    GGGGGG
                    GGGGGG
                    GGGGGG
                    HHHHHH
                    HHHHHH
                    HHHHHH"
        
        # Gather plots
        (model_out <- wrap_elements(study_plot+zoom_map+
                                      wrap_elements(add_ce_plot_ghm+add_ce_plot_sg)+ggtitle("Conditional Effects")+
                                      pp_dens+pp_err+
                                      trace+
                                      wrap_elements(full = param_table) + ggtitle("Model Coeeficient Table")+
                                      wrap_elements(full = data_sum_tbl) + ggtitle("Species Data Summary")+
                                      plot_layout(design = design))+
           ggtitle(glue("{sp} - {mod}")) + 
           theme(plot.title = element_text(size = 40, face = "bold")))
        
        # Write out plot
        ggsave(model_out, filename = glue("out/model_diagnostics/niche/{sp}.pdf"), 
               width = 10, height = 20, device = cairo_pdf)
        
      } # if add model is not NULL
    } else { #if the interaction CIs DONT overlap 0...
      
      #-- INTERACTIVE --#
      #- Model Basics -#        
      load(int_modlist_full[i]) # load model
      intmod <- out$model
      sp <- out$species
      out_int <- out
      fe_int <- fixef(out_int$model) #get fixed effects
      mod <- "Additive Model"
      
      data_sum_tbl <- out_int$data %>% 
        group_by(studyid) %>% 
        summarize(num_inds = n_distinct(ind_f),
                  #num_weeks = n_distinct(wk),
                  num_obs = n()) %>%
        mutate(obs_per_ind = num_obs/num_inds) %>% 
        adorn_totals() %>% 
        tableGrob()
      
      
      #- Make conditional predictions -#
      # get observed quantiles of ghm to set "low" and "high" human mod
      ghmq <- quantile(out$data$ghm_scale, probs = c(0.1, 0.9), na.rm = T)
      sgq <- quantile(out$data$sg_norm, probs = c(0.1, 0.9), na.rm = T)
      
      # Conditional Effects Plot for safegraph
      ce_int_sg <- conditional_effects(x=intmod, 
                                       effects = c("sg_norm:ghm_scale"),
                                       int_conditions = list(ghm_scale = ghmq),
                                       re_formula = NA, 
                                       plot = F) 
      # Plot
      (int_ce_plot_sg <-  plot(ce_int_sg, 
                               plot = F,
                               rug = F,
                               line_args = list("se" = T))[[1]] + 
          scale_color_manual(values = palnew, name = "Human \n Modification",
                             labels = c("High", "Low")) +
          scale_fill_manual(values = palnew, name = "Human \n Modification",
                            labels = c("High", "Low")) +
          xlab("Human Mobility") +
          ylab("Scaled Niche Breadth")+
          theme_minimal() +
          theme(panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.minor.x = element_blank(),
                # legend.position = "none",
                axis.title = element_text(size = 11,
                                          face = "bold"),
                axis.ticks = element_line(color = "#4a4e4d"),
                axis.text = element_text(size = 12),
                text = element_text(family = "Arial", color = "#4a4e4d")) 
      )
      
      # Conditional Effects Plot for GHM
      ce_int_ghm <- conditional_effects(x=intmod, 
                                        effects = c("ghm_scale:sg_norm"),
                                        int_conditions = list(sg_norm = sgq),
                                        re_formula = NA, 
                                        plot = F) 
      # Plot
      (int_ce_plot_ghm <-  plot(ce_int_ghm, 
                                plot = F,
                                rug = F,
                                line_args = list("se" = T))[[1]] + 
          scale_color_manual(values = palnew, name = "Human \n Mobility",
                             labels = c("High", "Low")) +
          scale_fill_manual(values = palnew, name = "Human \n Mobility",
                            labels = c("High", "Low")) +
          xlab("Human Modification") +
          ylab("Scaled Niche Breadth")+
          theme_minimal() +
          theme(panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.minor.x = element_blank(),
                # legend.position = "none",
                axis.title = element_text(size = 11,
                                          face = "bold"),
                axis.ticks = element_line(color = "#4a4e4d"),
                axis.text = element_text(size = 12),
                text = element_text(family = "Arial", color = "#4a4e4d")) 
      )
      
      #- Extract parameter table as grob -#
      param_table <- tableGrob(parameters(intmod))
      
      #- Posterior Predictive Plot -#
      (pp_dens <- pp_check(intmod)+ggtitle("Posterior Predictive Distribution"))
      
      #- Predictive Error Plot -#
      (pp_err <- pp_check(intmod, type='error_scatter_avg')+ggtitle("Posterior Predictive Errors"))  
      
      #- MCMC Trace Plot -#
      (trace <- mcmc_plot(intmod, type = "trace") +ggtitle("MCMC Traces"))    
      
      
      
      #-- Maps --#
      
      #- Daily Fixes Map -#
      
      # Filter to sp
      evt_sp <- evt0  %>% 
        filter(species==sp) %>% 
        mutate(doy = day(timestamp)) 
      
      # Summarize to daily fix, covert to sf
      evt_daily <- evt_sp %>% 
        st_as_sf(coords = c("lon", "lat"), crs=4326) %>% 
        group_by(study_id, ind_f, doy) %>% 
        summarize(st_union(geometry),
                  study_id = study_id[1]) %>% 
        st_centroid()
      
      # Define bounding range
      xrange <- c(min(evt_sp$lon)-1, max(evt_sp$lon)+1)
      yrange <- c(min(evt_sp$lat)-1, max(evt_sp$lat)+1)
      
      # Plot map
      (zoom_map <- ggplot()+
          geom_sf(data = us)+
          geom_sf(data = evt_daily, inherit.aes = F, aes(color = ind_f))+
          # geom_sf_label_repel(data = std_sf, aes(label = species), force_pull = 0, force = 10) +
          coord_sf() + 
          xlim(xrange) +
          ylim(yrange)+
          ggtitle("Individual Daily Positions")+
          theme(legend.position = "none") )
      
      #- Study Map -#
      
      # Get study centroids (based on points)
      std_sf <- evt_sp %>% 
        left_join(std0) %>% 
        st_as_sf(coords = c("lon", "lat"), crs=4326) %>% 
        group_by(study_id) %>% 
        summarize(st_union(geometry),
                  study_id = study_id[1],
                  study_name = study_name[1]) %>% 
        st_centroid()
      
      # Get bouding box of centroids to set plot lims
      stdbb <- st_bbox(std_sf)
      
      # make plot
      study_plot <- ggplot()+
        geom_sf(data = us)+
        geom_sf(data = std_sf, inherit.aes = F)+
        xlim(c(stdbb['xmin']-1, stdbb['xmax']+1)) +
        ylim(c(stdbb['ymin']-1, stdbb['ymax']+1))+
        ggtitle("Study Centroids") +
        geom_sf_label_repel(data = std_sf, aes(label = study_name), force_pull = 0, force = 10) +
        coord_sf()
      
      
      #---- Assemble Plots ---#
      
      # Design layout
      design <- "   
                  #IIII#
                  #IIII#
                  AAABBB
                  AAABBB
                  AAABBB
                  CCCCCC
                  CCCCCC
                  EEEFFF
                  EEEFFF
                  GGGGGG
                  GGGGGG
                  GGGGGG
                  GGGGGG
                  HHHHHH
                  HHHHHH
                  HHHHHH"
      
      # Gather plots
      (model_out <- wrap_elements(study_plot+zoom_map+
                                    wrap_elements(int_ce_plot_ghm+int_ce_plot_sg)+ggtitle("Conditional Effects")+
                                    pp_dens+pp_err+
                                    trace+
                                    wrap_elements(full = param_table) + ggtitle("Model Coeeficient Table")+
                                    wrap_elements(full = data_sum_tbl) + ggtitle("Species Data Summary")+
                                    plot_layout(design = design))+
         ggtitle(glue("{sp} - {mod}")) + 
         theme(plot.title = element_text(size = 40, face = "bold")))
      
      # Write out plot
      ggsave(model_out, filename = glue("out/model_diagnostics/niche/{sp}.pdf"), 
             width = 10, height = 20, device = cairo_pdf)
      
      
      
      
    } # else collect the interactions
  } else {#if int is NULL...
    #...then load the additive model instead
    if(add_modlist_full[i] != "NULL"){
      
      #- Model Basics -#        
      load(add_modlist_full[i]) # load model
      addmod <- out$model
      sp <- out$species
      out_add <- out
      fe_add <- fixef(out_add$model) #get fixed effects
      mod <- "Additive Model"
      
      data_sum_tbl <- out_add$data %>% 
        group_by(studyid) %>% 
        summarize(num_inds = n_distinct(ind_f),
                  #num_weeks = n_distinct(wk),
                  num_obs = n()) %>%
        mutate(obs_per_ind = num_obs/num_inds) %>% 
        adorn_totals() %>% 
        tableGrob()
      
      #- Make conditional predictions -#
      # get observed quantiles of ghm to set "low" and "high" human mod
      ghmq <- quantile(out$data$ghm_scale, probs = c(0.1, 0.9), na.rm = T)
      sgq <- quantile(out$data$sg_norm, probs = c(0.1, 0.9), na.rm = T)
      
      # Conditional Effects Plot for safegraph
      ce_add_sg <- conditional_effects(x=addmod, 
                                       effects = c("sg_norm"),
                                       int_conditions = list(ghm_scale = ghmq,
                                                             sg_norm = sgq),
                                       re_formula = NA, 
                                       plot = F) 
      # Plot
      (add_ce_plot_sg <-  plot(ce_add_sg, 
                               plot = F,
                               rug = F,
                               line_args = list("se" = T))[[1]] + 
          xlab("Human Mobility") +
          ylab("Scaled Niche Breadth")+
          theme_minimal() +
          theme(panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.minor.x = element_blank(),
                # legend.position = "none",
                axis.title = element_text(size = 11,
                                          face = "bold"),
                axis.ticks = element_line(color = "#4a4e4d"),
                axis.text = element_text(size = 12),
                text = element_text(family = "Arial", color = "#4a4e4d")) 
      )
      
      # Conditional Effects Plot for GHM
      ce_add_ghm <- conditional_effects(x=addmod, 
                                        effects = c("ghm_scale"),
                                        int_conditions = list(ghm_scale = ghmq,
                                                              sg_norm = sgq),
                                        re_formula = NA, 
                                        plot = F) 
      # Plot
      (add_ce_plot_ghm <-  plot(ce_add_ghm, 
                                plot = F,
                                rug = F,
                                line_args = list("se" = T))[[1]] + 
          xlab("Human Modification") +
          ylab("Scaled Niche Breadth")+
          theme_minimal() +
          theme(panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.minor.x = element_blank(),
                # legend.position = "none",
                axis.title = element_text(size = 11,
                                          face = "bold"),
                axis.ticks = element_line(color = "#4a4e4d"),
                axis.text = element_text(size = 12),
                text = element_text(family = "Arial", color = "#4a4e4d")) 
      )
      
      #- Extract parameter table as grob -#
      param_table <- tableGrob(parameters(addmod))
      
      #- Posterior Predictive Plot -#
      (pp_dens <- pp_check(addmod)+ggtitle("Posterior Predictive Distribution"))
      
      #- Predictive Error Plot -#
      (pp_err <- pp_check(addmod, type='error_scatter_avg')+ggtitle("Posterior Predictive Errors"))  
      
      #- MCMC Trace Plot -#
      (trace <- mcmc_plot(addmod, type = "trace") +ggtitle("MCMC Traces"))    
      
      
      
      #-- Maps --#
      
      #- Daily Fixes Map -#
      
      # Filter to sp
      evt_sp <- evt0  %>% 
        filter(species==sp) %>% 
        mutate(doy = day(timestamp)) 
      
      # Summarize to daily fix, covert to sf
      evt_daily <- evt_sp %>% 
        st_as_sf(coords = c("lon", "lat"), crs=4326) %>% 
        group_by(study_id, ind_f, doy) %>% 
        summarize(st_union(geometry),
                  study_id = study_id[1]) %>% 
        st_centroid()
      
      # Define bounding range
      xrange <- c(min(evt_sp$lon)-1, max(evt_sp$lon)+1)
      yrange <- c(min(evt_sp$lat)-1, max(evt_sp$lat)+1)
      
      # Plot map
      (zoom_map <- ggplot()+
          geom_sf(data = us)+
          geom_sf(data = evt_daily, inherit.aes = F, aes(color = ind_f))+
          # geom_sf_label_repel(data = std_sf, aes(label = species), force_pull = 0, force = 10) +
          coord_sf() + 
          xlim(xrange) +
          ylim(yrange)+
          ggtitle("Individual Daily Positions")+
          theme(legend.position = "none") )
      
      #- Study Map -#
      
      # Get study centroids (based on points)
      std_sf <- evt_sp %>% 
        left_join(std0) %>% 
        st_as_sf(coords = c("lon", "lat"), crs=4326) %>% 
        group_by(study_id) %>% 
        summarize(st_union(geometry),
                  study_id = study_id[1],
                  study_name = study_name[1]) %>% 
        st_centroid()
      
      # Get bouding box of centroids to set plot lims
      stdbb <- st_bbox(std_sf)
      
      # make plot
      study_plot <- ggplot()+
        geom_sf(data = us)+
        geom_sf(data = std_sf, inherit.aes = F)+
        xlim(c(stdbb['xmin']-1, stdbb['xmax']+1)) +
        ylim(c(stdbb['ymin']-1, stdbb['ymax']+1))+
        ggtitle("Study Centroids") +
        geom_sf_label_repel(data = std_sf, aes(label = study_name), force_pull = 0, force = 10) +
        coord_sf()
      
      
      #---- Assemble Plots ---#
      
      # Design layout
      design <- "   #GGGG#
                    #GGGG#
                    #AAAA#
                    #AAAA#
                    #AAAA#
                    BBBBBB
                    BBBBBB
                    CCCDDD
                    CCCDDD
                    EEEEEE
                    EEEEEE
                    EEEEEE
                    EEEEEE
                    FFFFFF
                    FFFFFF
                    FFFFFF"
      
      # Gather plots
      (model_out <- wrap_elements(study_plot + #A
                                    wrap_elements(add_ce_plot_ghm+add_ce_plot_sg)+ggtitle("Conditional Effects")+ #B
                                    pp_dens+pp_err+ #C, D
                                    trace+ #E
                                    wrap_elements(full = param_table) + ggtitle("Model Coeeficient Table")+ #F
                                    wrap_elements(full = data_sum_tbl) + ggtitle("Species Data Summary")+ #G
                                    plot_layout(design = design))+
         ggtitle(glue("{sp} - {mod}")) + 
         theme(plot.title = element_text(size = 40, face = "bold")))
      
      # Write out plot
      ggsave(model_out, filename = glue("out/model_diagnostics/niche/{sp}.pdf"), 
             width = 10, height = 20, device = cairo_pdf)
      
    } #fi
  } #else
  
}# i 

dbDisconnect(db)
#beepr::beep()
message("all done....")
