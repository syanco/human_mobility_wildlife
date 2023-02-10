library(dplyr)
library(brms)
library(ggplot2)

b0 <- 4
b1 <- 1.5
v <- 4
E <- 2

dat <- data.frame(x = rnorm(1000, b0, v)) %>% 
  mutate(y = (x*b1)+rnorm(1000, 0, E))

ggplot(dat, aes(x=x, y=y))+
  geom_point()

mod <- brm(
  y ~ x,
  data = dat,
  # family = poisson,
  init = 0,
  cores = 1,
  iter = 1000,
  thin = 1
)

mod