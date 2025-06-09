# Check the area and niche intra-individual interactive models autocorrelation
# This script runs the same checks twice, first for the area interactive model, 
# then for niche interactive model, and saves one plot per model 

library(tidyverse)

dir <- "~/repositories/human_mobility_wildlife/out/intra_ind_models/"

# --- AREA ---

load(file.path(dir, "area_interactive_rs-SY.rdata"))
test_dat <- test_mod$data

# Get posterior fitted values (mean)
test_dat$resid <- residuals(test_mod, summary = TRUE)[, "Estimate"]
test_dat$fitted <- fitted(test_mod, summary = TRUE)[, "Estimate"]

# Compute autocorrelation per individual
acf_plots <- test_dat %>%
  group_by(ind_f) %>%
  summarise(acf1 = acf(resid, 
                       plot = FALSE, 
                       lag.max = 1)$acf[2])  # lag-1 only

# Summarize across individuals
acf_summary <- acf_plots %>%
  summarise(mean_acf1 = mean(acf1, na.rm = TRUE),
            sd_acf1 = sd(acf1, na.rm = TRUE),
            n = n())

print(acf_summary)

plot <- ggplot(acf_plots, aes(x = acf1)) +
  geom_histogram(bins = 30, 
                 fill = "#3182bd", 
                 color = "white") +
  theme_minimal() +
  labs(x = "Lag-1 ACF of residuals", 
       y = "Count",
       title = "Area Intra-Ind Model: Residual autocorrelation across individuals")

ggsave(file.path(dir, "area_ac_check.png"))




# --- NICHE ---

load(file.path(dir, "niche_interactive_rs-SY.rdata"))
test_dat <- test_mod$data

# Get posterior fitted values (mean)
test_dat$resid <- residuals(test_mod, summary = TRUE)[, "Estimate"]
test_dat$fitted <- fitted(test_mod, summary = TRUE)[, "Estimate"]

# Compute autocorrelation per individual
acf_plots <- test_dat %>%
  group_by(ind_f) %>%
  summarise(acf1 = acf(resid, 
                       plot = FALSE, 
                       lag.max = 1)$acf[2])  # lag-1 only

# Summarize across individuals
acf_summary <- acf_plots %>%
  summarise(mean_acf1 = mean(acf1, na.rm = TRUE),
            sd_acf1 = sd(acf1, na.rm = TRUE),
            n = n())

print(acf_summary)

plot <- ggplot(acf_plots, aes(x = acf1)) +
  geom_histogram(bins = 30, 
                 fill = "#3182bd", 
                 color = "white") +
  theme_minimal() +
  labs(x = "Lag-1 ACF of residuals", 
       y = "Count",
       title = "Niche Intra-Ind Model: Residual autocorrelation across individuals")

ggsave(file.path(dir, "niche_ac_check.png"))
