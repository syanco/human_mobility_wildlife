library(brms)
library(tidyverse)

size <- read_csv("out/dbbmm_size.csv") %>% 
  mutate(trt_new = gsub('_.*','',trt),
         year_f = factor(year),
         trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld"))

form1 <- bf(area ~ trt_new*year_f + (1|species),
            ndt ~ trt_new*year_f + (1|species))
mod1 <- brm(form1, data = size, 
            family = shifted_lognormal(), cores = 2)
mod1
conditional_effects(mod1, points = T)
pp_check(mod1, type = "stat", stat = "mean")

mod2 <- brm(area ~ trt_new*year_f + (1+year_f|species), data = size, 
            family = shifted_lognormal(), cores = 2)
mod2
conditional_effects(mod2, points = T)
pp_check(mod2, type = "stat", stat = "mean")
pp_check(mod2) + xlim(c(0,100000))

mod3 <- brm(area ~ trt_new*year_f + (1+year_f*trt_new|species), data = size,
            family = lognormal(), cores = 2)
mod3
conditional_effects(mod3, points = T)
pp_check(mod3, type = "stat", stat = "mean")
pp_check(mod3) + xlim(c(0,100000))


mod1
conditional_effects(mod1, points = T)
pp_check(mod1, type = "stat", stat = "mean")

pp_check(mod1, type = "stat", stat = "mean") +
  coord_cartesian(xlim = c(0.001, 300000)) +
  scale_x_continuous("Area",
                     trans = "log",
                     breaks = c(0.001, 1, 100, 1000, 10000, 100000),
                     labels = c(
                       "0.001", "1", "100", "1000", "10000",
                       "100000"
                     )
  ) +
  ggtitle("Prior predictive distribution of means")

n <- size %>% 
  group_by(species) %>% 
  summarize(n = n())

