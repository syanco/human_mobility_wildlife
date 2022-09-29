#!/usr/bin/env Rscript 
#
#CONDA: covid
#
# This script generates individual dynamic brownian bridge models and associated 
# UDs for migratory  periods.

# TODO:  The dBBMM paramaters (e.g., window size, margin, error, etc.) are 
# currently hardcoded.  Could be passed in as options to the script.
# TODO: verify the volume/probability problem for write out.
# 
# ==== Setup ====

'
Estimate Dynamic Brownian Bridge Movement Models (dBBMs) for pre-segemented 
migratory periods across multiple individuals

Usage:
make_dbbmm.r <db> <out> <nc>
make_dbbmm.r (-h | --help)

Parameters:
  db: path to movement databse. 
  out: path to output directory.
  nc: number of cores for parallel processing
  
Options:
-h --help     Show this screen.
-v --version     Show version.
' -> doc

#---- Input Parameters ----#
if(interactive()) {
  library(here)
  
  .wd <- '~/projects/covid-19_movement'
  rd <- here::here
  
  .outPF <- file.path(.wd,'out')
  .dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
  
  .nc <- 2
  
} else {
  library(docopt)
  library(rprojroot)
  
  ag <- docopt(doc, version = '0.1\n')
  .wd <- getwd()
  
  source(file.path(.wd, 'analysis/src/funs/input_parse.r'))
  
  .outPF <- makePath(ag$out)
  .dbPF <- makePath(ag$db)
  .nc <- ag$nc
  
}

#---- Initialize Environment ----#
# Start time
t0 <- Sys.time()

# Run startup
source(file.path(.wd,'analysis/src/startup.r'))

# Load packages
suppressWarnings(
  suppressPackageStartupMessages({
    library(DBI)
    library(RSQLite)
    library(lubridate)
    # library(raster)
    # library(move)
    library(doMC)
    library(foreach)
    library(brms)
  }))
conflict_prefer("accumulate", "purrr")
conflict_prefer("ar", "brms")
conflict_prefer("lag", "stats")
conflict_prefer("when", "purrr")
#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

#---- Initialize database ----#
message("Initializing database connection...")

invisible(assert_that(file.exists(.dbPF)))
db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)
invisible(assert_that(length(dbListTables(db))>0))
indtb <- tbl(db,'individual') %>% 
  collect()
evttb <- tbl(db, 'event_clean') %>% 
  collect()

message("Disconnecting from databse...")
dbDisconnect(db)

# paired <- evttb %>% 
#   mutate(yrnum = as.numeric(yr)) %>% 
#   group_by(individual_id) %>% 
#   summarize(,
#             minyr = min(yrnum),
#             maxyr = max(yrnum),
#             paired = minyr != maxyr) %>% 
#   filter(paired)



evt_paired <- evttb %>% 
  filter(individual_id %in% paired_vec)

traits <- read_csv("raw_data/anthropause_data_sheet.csv")

size <- read_csv("out/dbbmm_size.csv") %>%
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>%
  distinct()

paired <- size %>% 
  # mutate(yrnum = as.numeric(yr)) %>% 
  group_by(ind_id) %>% 
  summarize(minyr = min(year),
            maxyr = max(year),
            paired = minyr != maxyr) %>% 
  filter(paired)

paired_vec <- paired %>% 
  pull(ind_id)

size_paired <- size %>% 
  filter(ind_id %in% paired_vec) %>% 
  left_join(traits, by = c("species" = "Species")) %>% 
  mutate(diet = case_when(Diet.PlantO >= 50 ~ "Herbivore",
                          Diet.Fruit >= 50 ~ "Frugivore",
                          Diet.Scav >= 50 ~ "Savenger",
                          Diet.Inv >= 50 ~ "Insectivore",
                          (Diet.Vend+Diet.Vect+Diet.Vfish) >= 50 ~ "Carnivore",
                          TRUE ~ "Omnivore"),
         log_area = log(area), #get log of weekly area use
         log_area_scale = scale(log_area), # standardize it
         sg_norm = sg / cbg_area, # normalize safegraph data by size of the CBG
         sg_sqrt = sqrt(sg_norm),
         sg_scale = scale(sg_norm),
         # log_sg_norm = log(sg_norm),
         ghm_scale = scale(ghm),
         ndvi_scale = scale(ndvi),
         lst_scale = scale(lst),
         ind_f = as.factor(ind_id), # create factor version of ind for REs
         grp = paste(ind_f, year, sep = "_"), # create indXyr grouping factor
         # trt_new = gsub('_.*','',trt),
         year_f = factor(year), # create year factor
         # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
         # sp2 = gsub(" ", "_", species),
         wk_n = as.numeric(substring(wk, 2)), # extract week number
         ts = parse_date_time(paste(year, wk, 01, sep = "-"), "%Y-%U-%u"), # make better date format
         study_f = as.factor(study_id),
         ind_wk = paste0(ind_id,wk))  # make study id factor) %>% 

size_mean <- size_paired %>% 
  group_by(ind_id, year_f) %>% 
  summarize(m_area = mean(area, na.rm = T))

size_wide <- size_paired %>% 
  pivot_wider(id_cols = c(ind_f, wk, species), 
              values_from = c(log_area_scale, sg_norm, ghm_scale), 
              names_from = year_f) %>% 
  mutate(area_diff = log_area_scale_2020-log_area_scale_2019,
         sg_diff = sg_norm_2020-sg_norm_2019,
         ghm_diff = ghm_scale_2020-ghm_scale_2019) %>% 
  filter(!is.nan(sg_diff)) %>% 
  filter(!is.na(sg_diff))

# get sample size per sp
size_wide %>% 
  group_by(species) %>% 
  summarize(nind = n_distinct(ind_id)) %>% 
  arrange(-nind)

#distributuion of delta sg
ggplot(size_wide)+
  geom_density(aes(x=sg_diff))+
  theme_minimal()+
  NULL
# change in SG across domain of GHM
ggplot(size_wide)+
  geom_point(aes(x=ghm_scale_2020, y = sg_diff))
  
ggplot(size_wide) +
  geom_line(aes(x=wk, y = diff, group = ind_id, color = log(sg_diff)), alpha = 0.3) +
  scale_color_viridis_c()+
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  facet_wrap(~species, scales = "free_y")


ggplot(size_wide,aes(x=wk, y = diff)) +
  geom_point(aes(color = log(sg_diff)), alpha = 0.3) +
  scale_color_viridis_c()+
  geom_smooth()+
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  facet_wrap(~species, scales = "free_y")

ggplot(size_wide,aes(x=sg_diff, y = area_diff)) +
  # geom_point(alpha=0.1)+
  geom_bin2d() +
  scale_fill_viridis_c() +
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  ylab("Change in Weekly Space Use \n (=2020-2019)")+
  xlab("Change in Human Mobility \n (=2020-2019)") +
  facet_wrap(~species, scales = "free")+
  theme_minimal()+
  theme(legend.position = "none") +
  NULL

ggplot(size_wide,aes(x=ghm_diff, y = area_diff)) +
  # geom_point()+
  geom_bin2d() +
  scale_fill_viridis_c() +
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  ylab("Change in Weekly Space Use \n (=2020-2019)")+
  xlab("Change in Human Modification \n (=2020-2019)") +
  facet_wrap(~species, scales = "free")+
  theme_minimal()+
  theme(legend.position = "none") +
  NULL



size_wide %>% 
  # filter(species == "Alces alces") %>% 
  ggplot(aes(x=sg_diff, y = area_diff)) +
  geom_point(aes(color = as.factor(ind_id)))+
  # geom_bin2d() +
  scale_fill_viridis_c() +
  geom_hline(aes(yintercept = 0), linetype = "dashed")+
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  ylab("Change in Weekly Space Use \n (=2020-2019)")+
  xlab("Change in Human Mobility \n (=2020-2019)") +
  # facet_wrap(~species, scales = "free")+
  theme_minimal()+
  theme(legend.position = "none") +
  NULL



# filter no diff approach


size_wide[size_wide$sg_diff == min(size_wide$sg_diff),]

x <- size_wide %>% 
  mutate(sg_diff2 = round(sg_diff, 5)) %>% 
  filter(sg_diff2 != 0)

t.test(size_mean$`2019`, size_mean$`2020`, paired = T)

ggplot(data = size_paired, aes(y = log_area_scale, x = year_f))+
  geom_boxplot()

ggplot(data = size_paired) +
  geom_line(aes(y = log_area_scale, x = year_f, group = ind_wk)) +
  facet_wrap(~species)

form <-  bf(area_diff ~ 1 + sg_diff2 + (1|species/ind_f) + ar(time = wk, gr = ind_f))
message("Fitting models with formula:")
print(form)

message("Starting model...")



# fit model
mod <- brm(
  form,
  data = x,
  family = student(),
  inits = 0,
  cores = 4,
  iter = 1000,
  thin = 2
)

mod
conditional_effects(mod)


ids <- unique(size_wide$ind_id)
v <- c()
for(i in 1:ids){
  v[i] <- size_wide %>% 
    filter(ind_id == ids[i]) %>% 
    summarize(v=var(sg_norm_2019)) %>% 
    pull(v)
}

(mv <- mean(v))
