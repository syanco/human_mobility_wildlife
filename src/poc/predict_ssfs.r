# scratch code to fit ssf and plot predictions...

library(brms)
library(dplyr)
library(INLA)
library(raster)

#loda data
dat <- read.csv("out/test_ssf_dat.csv")


mean.beta <- 0
prec.beta <- 1e-4  

formula.fixed <-  case ~  tmax_norm + 
  f(strt_n, model="iid", hyper = list(theta = list(initial = log(1e-6),
                                                    fixed=T))) +
  f(ind_f, tmax_norm, values=NULL, model="iid", 
    hyper=list(theta=list(initial=log(1), fixed=F,prior="pc.prec",
                          param=c(3,0.05))))

mod_inla <- inla(formula.fixed, family ="Poisson", data=dat,
                 control.fixed = list(
                   mean = mean.beta,
                   prec = list(default = prec.beta)))
summary(mod_inla)

# Predict on raster

r0 <- raster("raw_data/daymet_v4_tmax_annavg_hi_1991.tif")
plot(r0)

r_scale <- r0
values(r_scale) <- scale(values(r0))

nd_med <- median(values(r_scale), na.rm = T)
b1 <- mod_inla$summary.fixed["tmax_norm", "mean"]

rss_vec <- sapply(values(r_scale), FUN = function(x){
  (x*b1)-(nd_med*b1)
})

r_rss <- r_scale
values(r_rss) <- rss_vec
plot(r_rss)


########
tmax <- raster("raw_data/daymet_v4ll_daily_pr_tmax_202206.nc")
plot(tmax)


