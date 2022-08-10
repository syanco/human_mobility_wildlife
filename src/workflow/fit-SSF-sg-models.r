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
  .bgM <- file.path(.wd, "out/ssf-background-pts/moose")
  .outP <- file.path(.wd,'out/ssf-mods')
  
  .cores <- 20
  .minsp <- 10
  .iter <- 5000
  .thin <- 4
  
} else {
  library(docopt)
  
  ag <- docopt(doc, version = '0.1\n')
  
  .wd <- getwd()
  .cores <- ag$cores
  .minsp <- ag$minsp
  .iter  <- ag$iter
  .thin <- ag$thin
  
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
    library(amt)
    library(glmmTMB)
    library(RSQLite)
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
all_files <- list.files(.bgP)
# moose_id <- sub("\\..*", "", moose_files)
sum(str_detect(all_files, "ghm-sg"))
moose_ghm_files <- all_files[str_detect(all_files, "ghm-sg")]
str_extract(moose_ghm_files, "//d")


foreach(i = 1:nrow(sp_sum), .errorhandling = "pass", .inorder = F) %dopar% {
  
  # get focal species
  sp <- sp_sum$taxon_canonical_name[i]
  
  # list of individuals within species
  indls <- indtb %>% 
    filter(taxon_canonical_name == sp) %>% 
    pull(individual_id)
  
  # init empty list to stor annos
  sp_out <- list()
  
  # loop through individuals
  # j <- 2
  # 
  for(j in 1:length(moose_id)){
    ndvi <- read_csv(glue("{.bgP}/{moose_id[j]}_ndvi_1.csv")) %>% 
      mutate(lat = round(lat, 4),
             lng = round(lng, 4))      
    tmax <- read_csv(glue("{.bgP}/{moose_id[j]}_tmax_1.csv")) %>% 
      mutate(lat = round(lat, 4),
             lng = round(lng, 4))
    dat0 <- read_csv(glue("out/ssf-background-pts/individual-files/{moose_id[j]}.csv")) %>% 
      mutate(x2_ = round(x2_, 4),
             y2_ = round(y2_, 4))
    sg_ghm <- read_csv(glue("{.bgP}/individual-{moose_id[j]}ghm-sg.csv"))
    
    
    sp_out[[j]] <- dat0 %>% 
      left_join(ndvi, by = c("x2_" = "lng", "y2_" = "lat", 
                              "t2_" = "t2_",
                              "case_" = "case_")) %>% 
      left_join(tmax, by = c("x2_" = "lng", "y2_" = "lat", 
                             "t2_" = "t2_",
                             "case_" = "case_",
                             "crds" = "crds")) %>% 
      mutate(ind_id = indls[j],
             strt = paste0(ind_id,"_",burst_))
    
    # for(j in 1:20){
    #   # load the pts
    #   
    #   # Load data (and format for stoat)
    #   pts <- read_csv(glue("out/ssf-background-pts/individual-files/{indls[j]}.csv")) %>% 
    #     mutate(lng = x2_,
    #            lat = y2_,
    #            date = format(t2_, format = '%Y-%m-%d'))
    #   
    #   # annotate and store results
    #   message(glue("Annotating individual {j} of {length(indls)}"))
    #   sp_out[[j]] <- big_stoat(pts, vars) %>%
    #     select(-product) %>% 
    #     pivot_wider(names_from = variable, values_from = value, 
    #                 id_cols = c("x1_", "x2_", "y1_", "y2_", "sl_", "direction_p", 
    #                             "ta_", "t1_", "t2_", "dt_", "step_id_", "case_")) %>% 
    #     mutate(ind_id = indls[j])
  } # j
  
  
  dat0 <- do.call("rbind", sp_out)  
  
  dat <- dat0 %>%
    mutate(
      tmax_norm = scale(tmax),
      ndvi_norm = scale(ndvi)
    ) %>%
    distinct()    
  
  message(glue("Starting model for {sp}..."))
  
  # fit model
  TMBStruc.fix = glmmTMB(case_ ~ tmax_norm + ndvi_norm + (1|burst_) + 
                           (0 + tmax_norm | ind_id) + (0 + ndvi_norm | ind_id), 
                         family=poisson, data=dat, doFit=FALSE) 
  
  #' Then fix the standard deviation of the first random term
  TMBStruc.fix$parameters$theta[1] = log(1e3) 
  
  #' We need to tell `glmmTMB` not to change the variance by setting it to `NA`:
  TMBStruc.fix$mapArg = list(theta=factor(c(NA, 1:2))) # change the 2 to however many random slopes I have
  
  #' Then fit the model and look at the results:
  glmm.TMB.fixed = glmmTMB:::fitTMB(TMBStruc.fix) 
  summary(glmm.TMB.fixed)
  
  
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


