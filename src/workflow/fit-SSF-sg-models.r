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
  }))

# # Manage package conflicts
# conflict_prefer("accumulate", "foreach")
# conflict_prefer("ar", "brms")
# conflict_prefer("lag", "stats")
# conflict_prefer("when", "foreach")

# load breezy functions
source(file.path(.wd,'analysis/src/funs/auto/breezy_funs.r'))

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

# loop through species
foreach(i = 1:nrow(sp_sum), .errorhandling = "pass", .inorder = F) %dopar% {
  
  # get focal species
  sp <- sp_sum$taxon_canonical_name[i]
  
  
  #---- Initialize database ----#
  message(glue("Initializing database connection for individual {ind[j]}, year {yearvec[i]}..."))
  
  invisible(assert_that(file.exists(.dbPF)))
  db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
  invisible(assert_that(length(dbListTables(db))>0))
  
  message(glue("Starting model for {sp}..."))
  
  dat <- breadth %>%
    filter(scientificname == sp) %>% 
    mutate(
      sqrt_breadth = sqrt(total), #get log of weekly area use
      sqrt_breadth_scale = scale(sqrt_breadth), # standardize it
      sg_norm = scale(sg / cbg_area), # normalize safegraph data by size of the CBG
      # log_sg_norm = log(sg_norm),
      ghm_scale = scale(ghm),
      ndvi_scale = scale(ndvi.y),
      lst_scale = scale(lst.y),
      
      grp = paste(ind_f, year, sep = "_"), # create indXyr grouping factor
      # trt_new = gsub('_.*','',trt),
      year_f = factor(year), # create year factor
      # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
      # sp2 = gsub(" ", "_", species),
      wk_n = as.numeric(substring(week, 2)), # extract week number
      ts = parse_date_time(paste(year, week, 01, sep = "-"), "%Y-%U-%u"), # make better date format
      study_f = as.factor(studyid) # make study id factor
    ) %>%
    distinct()    
  
  # fit model
  mod <- brm(
    form,
    data = dat,
    # family = Gamma(link = "log"),
    inits = 0,
    cores = 1,
    iter = .iter,
    thin = .thin
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


