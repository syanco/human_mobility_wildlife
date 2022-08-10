library(ggplot2)
library(tidyverse)
# library(ggth)

moose_size <- read_csv("out/dbbmm_size.csv") %>% 
  # filter(species == "Haliaeetus leucocephalus") %>% 
  mutate(grp = paste0(ind_id, year))

moose_niche <- read_csv("out/niche_determinant_anthropause.csv") 

moose_size %>% 
  # filter(year == 2020) %>% 
  ggplot()+
  # geom_line(aes(x=wk, y = sg/cbg_area, color = as.factor(year), group = grp))+
  geom_smooth(aes(x=wk, y = sg/cbg_area)) +
  facet_wrap(~year)

moose_size %>% 
  # filter(year == 2020) %>% 
  ggplot()+
  geom_line(aes(x=wk, y = area, color = as.factor(year), group = grp))
  
  # geom_smooth(aes(x=wk, y = area)) +
  # facet_wrap(~year)

ggplot(moose_size)+
  geom_point(aes(x=sg/cbg_area, y = ghm, color = as.factor(year)))

ggplot(moose_size)+
  geom_density(aes(x=sg/cbg_area, group = as.factor(year), color = as.factor(year)))+
  xlim(c(0,0.0000008)) +
  theme_tufte()

length(unique(moose_niche$scientificname))
