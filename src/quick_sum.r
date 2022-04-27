library(ggplot2)
library(tidyverse)

load("out/quick_mod.rdata")


ggplot(size)+
  geom_density(aes(x=sg, color = year_f)) 
# +
#   xlim(c(0,5000000))

summary(size$area)

sum(size$area > 10000000)/nrow(size)

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


