# Check that area size estimates are not sample size contingent

library(tidyverse)
library(ggplot2)
library(lme4)
library(sjPlot)

size_dat <- read_csv("out/dbbmm_size.csv") %>% 
  # filter(area < 20000000000) %>% 
  mutate(ind_f = as.factor(ind_id),
         log_area = log(area), #get log of weekly area use
         # log_area_scale = scale(log_area), # standardize it
         area_scale = scale(area),
         sg_norm = scale(sg / cbg_area), # normalize safegraph data by size of the CBG
         # log_sg_norm = log(sg_norm),
         ghm_scale = scale(ghm),
         ndvi_scale = scale(ndvi),
         # lst_scale = scale(lst),
         tmax_scale = scale(tmax),
         # tmin_scale = scale(tmin),
         grp = paste(ind_f, year, sep = "_"), # create indXyr grouping factor
         # trt_new = gsub('_.*','',trt),
         year_f = factor(year), # create year factor
         # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
         # sp2 = gsub(" ", "_", species),
         wk_n = as.numeric(substring(wk, 2)), # extract week number
         ts = parse_date_time(paste(year, wk, 01, sep = "-"), "%Y-%U-%u"), # make better date format
         study_f = as.factor(study_id), # make study id factor
         # scale_n = scale(n)
  ) %>%
  distinct()  %>% 
  group_by(ind_f) %>%
  mutate(scale_n = scale(n),
         log_area_scale = scale(log_area)) %>%
  ungroup()
summary(size_dat$scale_n)


# ggplot(size_dat)+
#   geom_point(aes(x = n, y = area))+
#   facet_wrap(~species)
mod <- lmer(area ~ scale_n + (1|species/ind_f), data = size_dat)
mod_log<- lmer(log_area ~ scale_n + (1|species/ind_f), data = size_dat)

summary(mod)

confint(mod)

plot_model(mod, type = "eff")+
  theme_classic()+
  ylab("Area Size")+
  xlab("Sample Size (scaled)")+
  ggtitle("")



