library(ggplot2)
library(tidyverse)
library(ggthemes)
pal2 <- c("#E66100", "#5D3A9B")
dat <- read_csv("out/dbbmm_size.csv") %>% 
  # filter(species == sp[20]) %>%
  mutate(grp = paste0(ind_id, year))


ggplot(dat)+
  geom_density(aes(x=sg/cbg_area, group = as.factor(year), color = as.factor(year)),
               size = 2,show.legend=FALSE)+
  xlim(c(0,0.0000008)) +
  xlab("Human Mobility")+
  ylab("")+
  scale_color_manual(values = pal2, name = "Year", labels = c("2019", "2020")) +
  theme_tufte()+
  stat_density(aes(x=sg/cbg_area, colour=as.factor(year)),
               geom="line",position="identity", size = 0) +
  guides(colour = guide_legend(override.aes=list(size=1)))+
  theme(axis.line = element_line(size = 1)) + 
  facet_wrap(~species, scales = "free_y") +
  NULL
