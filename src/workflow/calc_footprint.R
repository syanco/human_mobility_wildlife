################################
####  Calc Area Footprint   ####
################################

# DESCRIPTION:  Calculates the "footprint of the anthropocene"


####---- INIT ----####

# Libraries
library(tidyverse)
library(glue)
library(wesanderson)

# Color Pallete
pal <- c("#EE7674", "#F4D5A4")


####---- Load Data ----####

area_pred_dat <- read_csv("out/area_change_prediction_2023-10-18.csv")
niche_pred_dat <- read_csv("out/niche_change_prediction_2023-10-19.csv")

birds <- c("Anas acuta", "Anas americana", "Anas clypeata", "Anas crecca",
           "Anas cyanoptera", "Anas platyrhynchos", "Anas strepera", 
           "Anser albifrons", "Anser caerulescens", "Aquila chrysaetos", 
           "Ardea alba", "Chen rossii", "Circus cyaneus", "Corvus corax", 
           "Haliaeetus leucocephalus", "Rallus longirostris")


####---- Calc Diffs ----####

# species list vector to loop over
area_spl <- unique(area_pred_dat$species)
niche_spl <- unique(niche_pred_dat$species)


##-- AREA --##
# init output list
area_diff_out <- list()

# loop over spp for area
for(i in 1:length(area_spl)){
  
  # Extract species data
  sp_dat <- area_pred_dat %>% 
    filter(species == area_spl[i])
  
  # Extract low estimate
  est_low <- sp_dat %>% 
    filter(ghm_case == "low" & sg_case == "low") %>% 
    pull(est_unscaled_exp)
  
  # Extract high estimate
  est_high <- sp_dat %>% 
    filter(ghm_case == "high" & sg_case == "high") %>% 
    pull(est_unscaled_exp)
  
  # Calc difference
  diff <- est_low-est_high
  
  # Stash outputs
  tmp_out <- tibble(species = area_spl[i],            # Species
                    est_low = est_low,           # Low Estimate
                    est_high = est_high,         # High Estimate
                    diff = diff,                 # Predicted difference
                    model = sp_dat$model[1],     # Model
                    tot_sig = sp_dat$tot_sig[1]) # Model significance flag
  
  area_diff_out[[i]] <- tmp_out 
}


#- Compile Results -#

# Convert results to df
area_diff_df <- do.call("rbind", area_diff_out) %>% 
  # filter(species != "Haliaeetus leucocephalus" ) %>% 
  mutate(diff_km = -diff/1000000, # change units
         prop = est_high/est_low, # get proportiona change
         perc_num = -round((1-prop)*100, 0), # percent change as a whole number
         perc_lab = glue("{perc_num}%")) %>% # percent change as text
  filter(tot_sig == "sig") %>% # only retain significant results
  mutate(direct = case_when(diff_km > 0 ~ "p", # estimate direction flag
                            diff_km < 0 ~ "n"))

# Calculate mean differnce
(area_sum_diff <- area_diff_df %>% 
    mutate(class = case_when(species %in% birds ~ "bird",
                             TRUE ~ "mammal")) %>% 
    group_by(class) %>% 
    summarize(meandiff = mean(diff_km),
              mindiff = min(diff_km),
              maxdiff = max(diff_km),
              meddiff = median(diff_km),
              meanprop = mean(prop),
              meanperc = mean(perc_num),
              medprop = median(prop),
              medperc = median(perc_num),
              maxperc = max(perc_num),
              minperc = min(perc_num)))

(area_grandsum_diff <- area_diff_df %>% 
    summarize(meandiff = mean(diff_km),
              mindiff = min(diff_km),
              maxdiff = max(diff_km),
              meddiff = median(diff_km),
              meanprop = mean(prop),
              meanperc = mean(perc_num),
              medprop = median(prop),
              medperc = median(perc_num),
              maxperc = max(perc_num),
              minperc = min(perc_num)))



##-- NICHE --##
# init output list
niche_diff_out <- list()

# loop over spp for area
for(i in 1:length(niche_spl)){
  
  # Extract species data
  sp_dat <- niche_pred_dat %>% 
    filter(species == niche_spl[i])
  
  # Extract low estimate
  est_low <- sp_dat %>% 
    filter(ghm_case == "low" & sg_case == "low") %>% 
    pull(est_unscaled_exp)
  
  # Extract high estimate
  est_high <- sp_dat %>% 
    filter(ghm_case == "high" & sg_case == "high") %>% 
    pull(est_unscaled_exp)
  
  # Calc difference
  diff <- est_low-est_high
  
  # Stash outputs
  tmp_out <- tibble(species = niche_spl[i],            # Species
                    est_low = est_low,           # Low Estimate
                    est_high = est_high,         # High Estimate
                    diff = diff,                 # Predicted difference
                    model = sp_dat$model[1],     # Model
                    tot_sig = sp_dat$tot_sig[1]) # Model significance flag
  
  niche_diff_out[[i]] <- tmp_out 
}


#- Compile Results -#

# Convert results to df
niche_diff_df <- do.call("rbind", niche_diff_out) %>% 
  # filter(species != "Haliaeetus leucocephalus" ) %>% 
  mutate(prop = est_high/est_low, # get proportional change
         perc_num = -round((1-prop)*100, 0), # percent change as a whole number
         perc_lab = glue("{perc_num}%")) %>% # percent change as text
  filter(tot_sig == "sig") %>% # only retain significant results
  mutate(direct = case_when(diff > 0 ~ "p", # estimate direction flag
                            diff < 0 ~ "n"))

# Calculate mean differnce
(niche_sum_diff <- niche_diff_df %>% 
    mutate(class = case_when(species %in% birds ~ "bird",
                             TRUE ~ "mammal")) %>% 
    group_by(class) %>% 
    summarize(meanprop = mean(prop),
              meanperc = mean(perc_num),
              medprop = median(prop),
              medperc = median(perc_num),
              maxperc = max(perc_num),
              minperc = min(perc_num)) )

(niche_grandsum_diff <- niche_diff_df %>% 
    summarize(meanprop = mean(prop),
              meanperc = mean(perc_num),
              medprop = median(prop),
              medperc = median(perc_num),
              maxperc = max(perc_num),
              minperc = min(perc_num)) )


# ####---- Plots ----####
# 
# # divergeing bar
# ggplot(diff_df) +
#   geom_bar(aes(y = reorder(species, diff_km), x = -diff_km, fill = direct, alpha = -perc_num), 
#            stat = "identity",
#            show.legend = FALSE) +
#   scale_fill_manual(values = rev(pal)) +
#   
#   geom_text(aes(y=species, x = -diff_km/2, label = perc_lab))+
#   geom_vline(aes(xintercept = -mean_diff), color = pal[1], linetype = "dashed")+
#   geom_vline(aes(xintercept = -0), color = "black")+
#   annotate(
#     geom = "text", 
#     x = -30, 
#     y = 9, 
#     label = glue("Mean Weekly Difference: {mean_diff} square km"),
#     vjust = 0.5,
#     hjust = 1,
#     color = "black"
#   ) +
#   annotate(
#     geom = "segment", 
#     x = -30, 
#     xend = -mean_diff, 
#     y = 9, 
#     yend = 9, 
#     arrow = arrow(length = unit(0.3, "cm")),
#     color = "black"
#   )+
#   ylab("Species") +
#   xlab(expression(paste(Delta, "Weekly Area Use (", km^2, ")"))) +
#   theme_minimal()+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())


