#!/usr/bin/env Rscript 
#
# ==== Input parameters ====

'
Usage:
fit-SSF-sg-models.r <dbPF> <bgP> <out> <cores> [<minsp> <iter> <thin>]
fit-SSF-sg-models.r (-h | --help)


Parameters:
nichedat: path to input niche breadth csv file.
vardat: path to dbbmm data csv with covars
out: path to output directory.
cores: the number of HPC cores to request

Options:
-h --help           Show this screen.
-v --version        Show version.

' -> doc

if(interactive()) {
  
  .wd <- getwd()
  
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
  .bgP <- file.path(.wd, "out/ssf-background-pts/annotated")
  # .bgM <- file.path(.wd, "out/ssf-background-pts/moose")
  .outP <- file.path(.wd,'out/ssf-mods')
  
  .cores <- 20
  .minsp <- 10
  .iter <- 5000
  .thin <- 4
  
} else {
  library(docopt)
  ag <- docopt(doc, version = '0.1\n')
  
  .wd <- getwd()
  
  source(file.path(.wd,'analysis/src/funs/input_parse.r'))
  
  .dbPF <- makePath(ag$dbPF)
  .bgP <- makePath(ag$bgP)
  # .bgM <- file.path(.wd, "out/ssf-background-pts/moose")
  .outP <- makePath(ag$out)
  
  .cores <- makePath(ag$cores)
  .minsp <- makePath(ag$minsp)
  .iter <- makePath(ag$iter)
  .thin <- makePath(ag$thin)
  
  source(file.path(.wd,'analysis/src/funs/input_parse.r'))
  
  .datPF <- makePath(ag$nichedat)
  .varPF <- makePath(ag$vardat)
  .outP <- makePath(ag$out)
  
}

# ==== Setup ====

#---- Initialize Environment ----#

t0 <- Sys.time()

source(file.path(.wd,'analysis/src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    # library(iterators)
    library(tidyverse)
    library(lubridate)
    library(glue)
    # library(foreach)
    # library(doMC)
    # library(amt)
    library(glmmTMB)
    library(RSQLite)
    library(broom)
    library(broom.mixed)
    # library(rstoat)
  }))

# # Manage package conflicts
# conflict_prefer("accumulate", "foreach")
# conflict_prefer("ar", "brms")
# conflict_prefer("lag", "stats")
# conflict_prefer("when", "foreach")

# load breezy functions
source(file.path(.wd,'analysis/src/funs/auto/breezy_funs.r'))
# source(file.path(.wd,'analysis/src/funs/big_stoat.r'))

# check arg inputs
.minsp <- ifelse(is.null(.minsp), 10, as.numeric(.minsp))
.iter <- ifelse(is.null(.iter), 3000, as.numeric(.iter))
.thin <- ifelse(is.null(.thin), 4, as.numeric(.thin))


#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))

indtb <- tbl(db,'individual') %>%
  collect()

message("Disconnecting from databse...")
dbDisconnect(db)

ind <- indtb %>%
  pull(individual_id)


#---- Load data ----#
# message("Loading data...")
# 
# dbbmms <- read_csv(.varPF)
# 
# # load and process data
# breadth <- read_csv(.datPF) %>%
#   filter(studyid != 351564596) %>%
#   filter(studyid != 1891587670) %>%
#   mutate(ind_f = as.factor(individual)) %>%  # create factor version of ind for REs)
#   left_join(dbbmms, by = c("scientificname" = "species", 
#                            "individual" = "ind_id", 
#                            "year" = "year", 
#                            "studyid" = "study_id", 
#                            "week" = "wk"))
# 
# get ind count per species
sp_sum <- indtb %>%
  group_by(taxon_canonical_name) %>%
  summarize(nind = length(unique(individual_id))) %>%
  filter(nind > .minsp) #require a minimum of 10 individuals

# ==== Start cluster and register backend ====
registerDoMC(.cores)

# ==== Perform analysis ====

# #declare model form
# form <-  bf(sqrt_breadth ~ sg_norm + ndvi_scale + lst_scale + (1 |grp) + ar(time = week, gr = grp))
# message("Fitting models with formula:")
# print(form)
# vars <- list(list(id = "MODIS/061/MOD11A1",
#                   static = FALSE,
#                   reducers = list("mean"),
#                   s_buff = 1000,
#                   t_buff = 1,
#                   bands = list("LST_Day_1km")),
#              list(id = "MODIS/MOD09GA_006_NDVI",
#                   static = FALSE,
#                   reducers = list("mean"),
#                   s_buff = 1000,
#                   t_buff = 1,
#                   bands = list("NDVI")))

# loop through species
# i <- 1
# 
# 
# all_files <- list.files(.bgP)
# # moose_id <- sub("\\..*", "", moose_files)
# sum(str_detect(all_files, "ghm-sg"))
# moose_ghm_files <- all_files[str_detect(all_files, "ghm-sg")]
# str_extract(moose_ghm_files, "//d")


foreach(i = 1:nrow(sp_sum), .errorhandling = "pass", .inorder = F) %dopar% {
  
  # get focal species
  sp <- sp_sum$taxon_canonical_name[i]
  
  # get list of ind for ssf
  fl <- list.files(glue("{.bgP}"))
  
  # list of individuals within species
  # indls <- indtb %>% 
  #   filter(taxon_canonical_name == sp) %>% 
  #   pull(individual_id)
  
  # init empty list to stor annos
  sp_out <- list()
  
  # loop through individuals
  # j <- 2
  # 
  for(j in 1:length(fl)){
    #get ind id
    ind <- str_extract(fl[j], "[^.]+")
    
    
    o <- tryCatch({
      ndvi <- read_csv(glue("{.bgP}/{ind}_ndvi_1.csv")) %>% 
        mutate(lat = round(lat, 4),
               lng = round(lng, 4))      
      tmax <- read_csv(glue("{.bgP}/{ind}_tmax_1.csv")) %>% 
        mutate(lat = round(lat, 4),
               lng = round(lng, 4))
      ghm_sg <- read_csv(glue("{.bgP}/moose/individual-{ind}-sg-ghm.csv"))%>% 
        mutate(x2_ = round(x2_, 4),
               y2_ = round(y2_, 4))
      dat0 <- read_csv(glue("{.bgP}/{ind}.csv")) %>% 
        mutate(x2_ = round(x2_, 4),
               y2_ = round(y2_, 4))
      
      
      sp_out[[j]] <- dat0 %>% 
        left_join(ndvi, by = c("x2_" = "lng", "y2_" = "lat", 
                               "t2_" = "t2_",
                               "case_" = "case_")) %>% 
        left_join(tmax, by = c("x2_" = "lng", "y2_" = "lat", 
                               "t2_" = "t2_",
                               "case_" = "case_",
                               "crds" = "crds")) %>%       
        left_join(ghm_sg, by = c("x2_" = "x2_", "y2_" = "y2_", 
                                 "t2_" = "t2_",
                                 "case_" = "case_",
                                 "burst_" = "burst_")) %>% 
        mutate(ind_id = ind,
               strt = paste0(ind,"_",step_id_.x))
    },error=function(e) e)
    
    
    if(inherits(o, "error")) next
    
  } # j
  
  
  dat0 <- do.call("rbind", sp_out)  
  
  dat <- dat0 %>%
    filter(ndvi > 0 & ndvi <1) %>% 
    mutate(
      tmax_norm = scale(tmax),
      ndvi_norm = scale(ndvi),
      cbg_km = cbg_area_m2/1000000,
      sg_norm = safegraph_daily_count/cbg_km,
      ghm_norm = scale(ghm),
      ind_f = as.numeric(as.factor(ind_id)),
      an_ID = as.numeric(ind_f),
      burst_f = as.factor(burst_),
      strt_n = as.numeric(as.factor(strt)),
      case = as.numeric(case_),
      yr = as.factor(year(t2_))
      
    ) %>%
    distinct()    
  
  
  message(glue("Starting model for {sp}..."))
  
  # fit SG model
  mod_sg = glmmTMB(case ~ sg_norm + 
                     ndvi_norm + tmax_norm +
                     (1|strt_n) +
                   (0 + sg_norm| ind_f) +
                   (0 + tmax_norm | ind_f) + (0 + ndvi_norm | ind_f),
                   family=poisson, data=dat, doFit=FALSE) 
  
  #' Then fix the standard deviation of the first random term
  mod_sg$parameters$theta[1] = log(1e3) 
  
  #' We need to tell `glmmTMB` not to change the variance by setting it to `NA`:
  mod_sg$mapArg = list(theta=factor(c(NA, 1:3))) # change the 2 to however many random slopes I have
  
  #' Then fit the model and look at the results:
  fit_mod_sg = glmmTMB:::fitTMB(mod_sg) 
  
  (sg_t <- broom.mixed::tidy(fit_mod_sg, effects = "fixed", conf.int = T,
                             exponentiate = T))
  # ggplot(sg_t) +
  #   geom_point(aes(y=term, x = estimate))+
  #   geom_errorbar(aes(y=term, xmin = conf.low, xmax= conf.high), width = 0.25)+
  #   geom_vline(aes(xintercept = 1), linetype = "dashed") +
  #   theme_minimal()
  
  # fit GHM model
  mod_ghm = glmmTMB(case_ ~ ghm_norm + 
                      tmax_norm + ndvi_norm +
                      (1|strt_n) + 
                      (0 + ghm_norm | ind_f)+
                      (0 + tmax_norm | ind_f) + (0 + ndvi_norm | ind_f),
                    family=poisson, data=dat, doFit=FALSE) 
  
  #' Then fix the standard deviation of the first random term
  mod_ghm$parameters$theta[1] = log(1e3) 
  
  #' We need to tell `glmmTMB` not to change the variance by setting it to `NA`:
  mod_ghm$mapArg = list(theta=factor(c(NA, 1:3))) # change the 2 to however many random slopes I have
  
  #' Then fit the model and look at the results:
  fit_mod_ghm = glmmTMB:::fitTMB(mod_ghm) 
  (ghm_t <- broom.mixed::tidy(fit_mod_ghm, effects = "fixed", conf.int = T,
                              exponentiate = T))
  
  #stash results into named list
  out <- list(
    species = sp,
    data = dat,
    sg_model = mod_sg,
    ghm_model = mod_ghm,
    sg_tidy = sg_t,
    ghm_tidy = ghm_t
  )
  
  #write out results
  save(out, file = glue("{.outP}/{sp}_SSF_out_{Sys.Date()}.rdata"))
  
} # i

# ==== Finalize script ====
message(glue('Script complete in {diffmin(t0)} minutes'))


