library(ggplot2)
library(tidyverse)
library(brms)
library(tidybayes)    # Manipulate Stan objects in a tidy way
library(broom)        # Convert model objects to data frames
library(broom.mixed)  # Convert brms model objects to data frames
library(emmeans)
library(lubridate)

load("out/quick_mod_rslopes.rdata")

mod <- mod_res_rslopes
mod
size <- read_csv("out/dbbmm_size.csv") %>% 
  mutate(
    # trt_new = gsub('_.*','',trt),
    year_f = factor(year),
    # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
  ) %>%
  mutate(sp2 = gsub(" ", "_", species)) %>% 
  # also remove any species not in tree (i.e. NAs)
  # filter(sp2 %in% rownames(phylo_vcov)) %>% 
  mutate(ind_f = as.factor(ind_id),
         wk_n = as.numeric(substring(wk,2)),
         ts = parse_date_time(paste(year, wk, 01, sep = "-"), "%Y-%U-%u" )) %>%
  distinct()



#- Load ain the species trait data
traits <- read_csv("raw_data/anthropause_data_sheet.csv") %>% 
  mutate(mig_mod = case_when(migratory == "Partial" ~ "migratory",
                             .$migratory == "non-migratory" ~ "non-migratory",
                             .$migratory == "Complete" ~ "migratory",
                             .$migratory == "migratory" ~ "migratory",
                             .$migratory == "Semi-nomadic" ~ "non-migratory",
                             .$migratory == "Unknown" ~ NA_character_))

size <- size %>% 
  left_join(traits, by = c("species" = "Species"))

size <- size %>%
  filter(mig_mod == "non-migratory") %>%
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>%
  mutate(log_area = log(area),
         sg_norm = area/cbg_area,
         log_sg_norm = log(sg_norm),
         grp = as.factor(paste(ind_f, year, sep = "_")),
         year_f = as.factor(year)) %>%
  group_by(grp) %>%
  arrange(.by_group = T)

# mod_res2 %>% 
#   emmeans(~ sg_norm,
#           at = list(ghm = 0.5),
#           epred = TRUE)
(trnd <- emtrends(mod, ~sg_norm*ghm, var = "sg_norm",
                  at = list(ghm = seq(0, 1, by = 0.1)),
                  epred = T, re_formula = NA) %>%
    gather_emmeans_draws())

(plot_grand_mean <- ggplot(trnd,
                           aes(x = .value / 10, fill = factor(ghm))) +
    stat_halfeye(slab_alpha = 0.75) +
    # scale_fill_okabe_ito(order = c(3, 4)) +
    ggtitle("Marginal effect of human mobility on weekly space use")+
    labs(x = "Average marginal effect of a\n0.1-point increase in the GHM",
         y = "Density", fill = "GHM") +
    # theme_clean() + 
    theme(legend.position = "bottom"))



ce <- conditional_effects(x=mod, 
                    effects = "sg_norm:ghm",
                    int_conditions = list(ghm = seq(0, 1, by = 0.25)),
                    # rug = T,
                    theme = theme_minimal(),
                    re_formula = NA)

ce_plot <-  plot(ce, plot = FALSE)[[1]] +
  theme_minimal()

ce_plot

emmip(mod_res2, ~sg_norm*ghm,
      at = list(ghm = seq(0, 1, by = 0.1)))



grand_mean_mod <- mod_res2 %>% 
  epred_draws(newdata = expand_grid(
    sg_norm = seq(0, 500, by = 10),
    ghm = seq(0, 1, by = 0.1), 
    # ndvi = seq(0, 1, by = 0.1),
    ndvi = mean(size$ndvi),
    # lst = seq(12400, 16250, by = 10),
    lst = mean(na.omit(size$lst)),
    grp = levels(size$grp),
    wk = seq(1, 35, by = 1)), 
    re_formula = NA)




nd <- expand_grid(
  sg_norm = seq(0, 500, by = 10),
  ghm = seq(0, 1, by = 0.1), 
  # ndvi = seq(0, 1, by = 0.1),
  ndvi = mean(size$ndvi),
  # lst = seq(12400, 16250, by = 10),
  lst = mean(na.omit(size$lst)),
  wk = seq(1, 35, by = 1),
  # wk = 1,
  grp = levels(size$grp),
  sp2 = unique(size$sp2))

x <- fitted(mod_res2, allow_new_levels = F) %>%
  data.frame()

size_res %>%
  data_grid(sg_norm,)

ggplot(data = size_res, aes(x = sg_norm)) +
  geom_smooth(data = x,
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5),
              stat = "identity",
              alpha = 1/4, size = 1/2)

pred <- mod_res2 %>% 
  epred_draws(newdata = nd,
              newdata2 = list(wk = rep(seq(1, 35, by = 1), 561)),
              re_formula = NULL,
              allow_new_levels = T)

plot_grand_mean_trend <- ggplot(grand_mean_civlib_dist, 
                                 aes(x = civil_liberties, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Civil liberties index", y = "Predicted media freedom index",
       fill = "Credible interval") +
  theme_clean() +
  theme(legend.position = "bottom")




size %>%
  add_epred_draws(mod_res2) %>%
  ggplot(aes(x = sg_norm, y = mpg, color = ordered(cyl))) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = mtcars) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")


ggplot(pred, aes(x = civil_liberties, y = .epred)) +
  stat_lineribbon() + 
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Civil liberties index", y = "Predicted media freedom index",
       fill = "Credible interval") +
  theme_clean() +
  theme(legend.position = "bottom")



################################################
ggplot(size_res)+
  geom_density(aes(x=sg_norm, color = year_f)) 
# +
#   xlim(c(0,5000000))

summary(size_res$area)

sum(size_res$area > 10000000)/nrow(size_res)

ggplot(size)+
  geom_point(aes(x=sg, y = area, color = ind_f)) +
  theme(legend.position = "none")

ggplot(size)+
  geom_point(aes(x=sg, y = pop, color = year_f)) 


size$ind_f[size$area == max(size$area)]

theme(legend.position = "none")


sum(size$area <1000)


############################
# Plots to diagnose dBBMMs #
############################

size <- read.csv("out/dbbmm_size.csv")

size %>% 
  filter(n > 1000) %>%
  summarise(max(area))


# 
# %>%
#   ggplot() +
#   geom_density(aes(x = area))



##-- Area ~ n(evts) --##

size %>% 
  # filter(mig_mod == "non-migratory") %>%
  # summarise(n())
  ggplot()+
  geom_point(aes(x=n, y = area))



maxsizesum <- size %>%
  group_by(species) %>%
  summarise(max(area),
            n(),
            mean(n)) %>%
  arrange(desc(`max(area)`))


size %>%
  filter(ind_id !=  size$ind_id[size$area == max(size$area)]) %>%
  ggplot()+
  geom_point(aes(x=n, y = area, color = as.factor(ind_id)))+
  xlim(0,500)+
  theme(legend.position = "none")


ggplot(size)+
  geom_point(aes(x=n, y = a_bb))


ggplot(size)+
  geom_point(aes(x=a_bb, y = area))

ggplot(size)+
  geom_point(aes(x=fixmed, y = area))

ggplot(size)+
  geom_point(aes(x=m_error, y = area))


size$ind_id[size$area == max(size$area)]

size %>%
  group_by(as.factor(ind_id))%>%
  summarize(m = max(area))%>%
  arrange(desc(m))

size%>% 
  # filter(n > 500)%>%
  summarise(n())

ggplot(size) +
  geom_point(aes(x=pop, y = sg))


tab <- size %>%
  group_by(ind_f) %>%
  summarize(mean_area = mean(area),
            sp = sp2[1]) %>%
  arrange(desc(mean_area))

write_csv(tab, file = "out/ind_area_04112022.csv")


