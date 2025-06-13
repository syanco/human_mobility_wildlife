# Check that individual weekly utilization distribution area size estimates are 
# not sample size contingent. Plot the relationship between scaled sample size 
# and area. Use output as supplementary figure.

library(tidyverse)
library(ggplot2)
library(lme4)
library(sjPlot)
library(here)

size_dat <- read_csv(here("out/dbbmm_size.csv")) %>% 
  # filter(area < 20000000000) %>% 
  mutate(ind_f = as.factor(ind_id),
         log_area = log(area), # get log of weekly area use
         area_scale = scale(area),
         sg_norm = scale(sg / cbg_area), # normalize safegraph data by CBG size 
         ghm_scale = scale(ghm),
         ndvi_scale = scale(ndvi),
         tmax_scale = scale(tmax),
         grp = paste(ind_f, year, sep = "_"), # create indXyr grouping factor
         year_f = factor(year), # create year factor
         wk_n = as.numeric(substring(wk, 2)), # extract week number
         ts = parse_date_time(paste(year, wk, 01, 
                                    sep = "-"), "%Y-%U-%u"), # date formatting
         study_f = as.factor(study_id),
         area_km2 = area/1000000) %>% # make study id factor
  distinct()  %>% 
  group_by(ind_f) %>%
  mutate(scale_n = scale(n),
         log_area_scale = scale(log_area)) %>%
  ungroup()

summary(size_dat$scale_n)

mod <- lmer(area_km2 ~ scale_n + (1|species/ind_f), data = size_dat)

summary(mod)

confint(mod)

plot_model(mod, 
           type = "eff",
           terms = "scale_n")+
  theme_classic()+
  ylab("Area Size (kilometers²)")+
  xlab("Sample Size (scaled)")+
  ggtitle("")

ggsave(here("out/check_area_sample_size_balance_dim9-6.png"), 
       width = 9, 
       height = 6)



