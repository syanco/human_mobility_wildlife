library(brms)
library(dplyr)
library(survival)
library(INLA)


# sim code from:  https://github.com/glmmTMB/glmmTMB/issues/535
group_effects <- rnorm(10,sd=2)
dat0 <- data.frame(group=rep(1:10,each=20*21),
                   var1=rnorm(4200),
                   stratum=rep(1:200,each=21)) %>%
  group_by(stratum) %>% 
  mutate(prob = exp(var1*.1+group_effects[group]),
         prob=prob/sum(prob))
dat0$used <- unlist(by(dat0,
                       factor(dat0$stratum),
                       function(x) rmultinom(1,1,prob=x$prob),
                       simplify=T))

# make new data for predict (we do this heer b/c inla i weird)
nd1 <- dat0 %>% 
  ungroup() %>% 
  summarize(stratum= 20,
            var1 = median(var1))
nd2 <- dat0 %>%
  ungroup() %>% 
  summarize(stratum = 1,
            var1 = quantile(var1, 0.05))

dat1 <- dat0 %>% filter(group == 1)

# fit the clogit model first
mod_cl <- clogit(used ~ var1 + strata(stratum), data=dat1) 


# fit INLA IPP
inla_dat <- data.frame(used = c(NA, NA),
                       stratum = c(NA,NA),
                       group = NA,
                       var1 = c(nd1$var1, nd2$var1))
dat_inla <- rbind(dat1, inla_dat)

mean.beta <- 0
prec.beta <- 1e-4  

formula.fixed <-  0 + used ~  var1 + 
  f(stratum, model="iid", hyper = list(theta = list(initial = log(1e-6),
                                                    fixed=T))) 

mod_inla <- inla(formula.fixed, family ="Poisson", data=dat_inla,
                 control.fixed = list(
                   mean = mean.beta,
                   prec = list(default = prec.beta)),
                 control.predictor = list(compute=TRUE, link = 1))


# # fit the IPP models next
# get_prior(used ~ -1 + var1 + (1|stratum), dat0)
# mod_ipp <- brm(
#   used ~ -1 + var1 + (1|stratum),
#   data = dat0,
#   family = poisson,
#   init = 0,
#   cores = 4,
#   iter = 1000,
#   thin = 1,
#   prior = set_prior("constant(1000000)", class = "sd", group = "stratum", coef = "Intercept")
# )

summary(mod_cl)
# summary(mod_ipp)
summary(mod_inla)

#CL
uncenter <- sum(coef(mod_cl) * mod_cl$means, na.rm=TRUE)
(p_1 <- predict(mod_cl, nd1, type = "lp", reference = "sample",
             se.fit = TRUE))
p1 <- p_1$fit+uncenter

(p_2 <- predict(mod_cl, nd2, type = "lp", reference = "sample",
              se.fit = TRUE))
p2 <- p_2$fit+uncenter

p1-p2

# (p <- predict(mod_ipp, nd1, re_formula = NA))

b <- mod_inla$summary.fixed["var1", "mean"]

unname(b*nd1$var1-b*nd2$var1)


