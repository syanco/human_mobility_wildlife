#!/usr/bin/env Rscript --vanilla

# ==== Input parameters ====

'
Usage:
fit-space-use-models.r <dat> <out> [--cores=<cores>] [--minsp=<minsp>] [--iter=<iter>]
fit-space-use-models.r (-h | --help)


Parameters:
dat: path to input csv file.
out: path to output directory.

Options:
-h --help           Show this screen.
-v --version        Show version.
-c --cores=<cores>  The number of cores
-m --minsp=<minsp>  Minimum number of indibiduals per species to run a model
-i --iter=<iter>    MCMC iterations
' -> doc

if(interactive()) {

  .wd <- getwd()
  
  .datPF <- file.path(.wd,'out/dbbmm_size.csv')
  .outP <- file.path(.wd,'out/single_species_models/area')

  .cores <- 1
  .minsp <- 10
  .iter <- 3000
  
} else {
  library(docopt)

  ag <- docopt(doc, version = '0.1\n')
  
  .wd <- getwd()
  .cores <- ag$cores
  .minsp <- ag$minsp
  .iter  <- ag$iter
  
  source(file.path(.wd,'analysis/src/funs/input_parse.r'))
  
  .datPF <- makePath(ag$dat)
  .outP <- makePath(ag$out)
  
}

# ==== Setup ====

#---- Initialize Environment ----#

t0 <- Sys.time()

source(file.path(.wd,'analysis/src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    # library(iterators)
    library(brms)
    library(tidyverse)
    library(lubridate)
    library(glue)
    library(foreach)
    library(doMC)
    
  }))

# Manage package conflicts
conflict_prefer("accumulate", "foreach")
conflict_prefer("ar", "brms")
conflict_prefer("lag", "stats")
conflict_prefer("when", "foreach")

source(file.path(.wd,'analysis/src/funs/auto/breezy_funs.r'))

.minsp <- ifelse(is.null(.minsp), 10, as.numeric(.minsp))
.iter <- ifelse(is.null(.iter), 3000, as.numeric(.iter))

#---- Local parameters ----#

#---- Files and directories ----#


#---- Load data ----#
message("Loading data...")

size <- read_csv("out/dbbmm_size.csv") %>%
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>%
  mutate(
    log_area = log(area),
    sg_norm = area / cbg_area,
    log_sg_norm = log(sg_norm),
    ind_f = as.factor(ind_id),
    grp = paste(ind_f, year, sep = "_"),
    # trt_new = gsub('_.*','',trt),
    year_f = factor(year),
    # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
    # sp2 = gsub(" ", "_", species),  
    wk_n = as.numeric(substring(wk, 2)),
    ts = parse_date_time(paste(year, wk, 01, sep = "-"), "%Y-%U-%u"),
    study_f = as.factor(study_id)
  ) %>%
  distinct()

# get ind count per species
sp_sum <- size %>%
  group_by(species) %>%
  summarize(nind = length(unique(ind_f))) %>%
  filter(nind > .minsp) #require a minimum of 10 individuals

# ==== Start cluster and register backend ====
registerDoMC(.cores)

# ==== Perform analysis ====

#declare model form
form <-  bf(log_area ~ sg_norm*ghm + ndvi + lst + (1 |study_f/grp) + ar(time = wk, gr = grp))
message("Fitting models with formula:")
print(form)

# loop throuigh species
foreach(i = 1:nrow(sp_sum), .errorhandling = "pass", .inorder = F) %dopar% {
  
  # get focal species
  sp <- sp_sum$species[i]
  
  message(glue("Strating model for {sp}..."))
  
  #filter data
  dat <- size %>%
    filter(species == sp)
  
  # fit model
  mod <- brm(
    form,
    # data2 = list(phylo_vcov = phylo_vcov),
    data = dat,
    # family = Gamma(link = "log"),
    inits = 0,
    cores = 1,
    iter = .iter,
    thin = 4
  )
  
  #stash results into named list
  out <- list(
    species = sp,
    data = dat,
    model = mod
  )
  
  #write out results
  save(out, file = glue("{.outP}/{sp}_{Sys.Date()}.rdata"))

} # i

# ==== Finalize script ====
message(glue('Script complete in {diffmin(t0)} minutes'))


