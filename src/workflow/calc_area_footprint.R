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

pred_dat <- read_csv("out/area_change_prediction_2023-06-23.csv")


####---- Calc Diffs ----####

# species list vector to loop over
spl <- unique(pred_dat$species)

# init output list
diff_out <- list()

# loop over spp
for(i in 1:length(spl)){
  
  # Extract species data
  sp_dat <- pred_dat %>% 
    filter(species == spl[i])
  
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
  tmp_out <- tibble(species = spl[i],            # Species
                    est_low = est_low,           # Low Estimate
                    est_high = est_high,         # High Estimate
                    diff = diff,                 # Predicted difference
                    model = sp_dat$model[1],     # Model
                    tot_sig = sp_dat$tot_sig[1]) # Model significance flag
  
  diff_out[[i]] <- tmp_out 
}


####---- Compile Results ----####

# Convert results to df
diff_df <- do.call("rbind", diff_out) %>% 
  mutate(diff_km = diff/1000000, # change units
         prop = est_high/est_low, # get proportiona change
         perc_num = -round((1-prop)*100, 0), # percent change as a whole number
         perc_lab = glue("{perc_num}%")) %>% # percent change as text
  filter(tot_sig == "sig") %>% # only retain significant results
  mutate(direct = case_when(diff_km > 0 ~ "p", # estimate direction flag
                            diff_km < 0 ~ "n"))

# Calculate mean differnec
mean_diff <- diff_df %>% 
  summarize(mdiff = mean(diff_km)) %>% 
  pull(mdiff) %>% 
  round(1)

med_diff <- diff_df %>% 
  summarize(mdiff = median(diff_km)) %>% 
  pull(mdiff) %>% 
  round(1)

####---- Plots ----####

# divergeing bar
ggplot(diff_df) +
  geom_bar(aes(y = reorder(species, diff_km), x = -diff_km, fill = direct, alpha = -perc_num), 
           stat = "identity",
           show.legend = FALSE) +
  scale_fill_manual(values = rev(pal)) +
  
  geom_text(aes(y=species, x = -diff_km/2, label = perc_lab))+
  geom_vline(aes(xintercept = -mean_diff), color = pal[1], linetype = "dashed")+
  geom_vline(aes(xintercept = -0), color = "black")+
  annotate(
    geom = "text", 
    x = -30, 
    y = 9, 
    label = glue("Mean Weekly Difference: {mean_diff} square km"),
    vjust = 0.5,
    hjust = 1,
    color = "black"
  ) +
  annotate(
    geom = "segment", 
    x = -30, 
    xend = -mean_diff, 
    y = 9, 
    yend = 9, 
    arrow = arrow(length = unit(0.3, "cm")),
    color = "black"
  )+
  ylab("Species") +
  xlab(expression(paste(Delta, "Weekly Area Use (", km^2, ")"))) +
  theme_minimal()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


