library(ggplot2)

ggplot(size)+
  geom_density(aes(x=area, color = mig_mod)) +
  xlim(c(0,5000000))

summary(size$area)

sum(size$area > 10000000)/nrow(size)
