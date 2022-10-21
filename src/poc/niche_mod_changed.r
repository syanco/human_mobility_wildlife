# What changed in the niche model???
# 
# 
# ++++++++++++++++++++++++++++++++++++++++++


.wd <- getwd()

.outPF <- file.path(.wd,'out/niche_determinant_anthropause.csv')
.dbPF <- file.path(.wd,'processed_data/mosey_mod.db')
.varPF <- file.path(.wd, "out/dbbmm_size.csv")
.iter <- 10000
.thin <- 5



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
    # library(doMC)
    # library(foreach)
    require(MVNH)
    library(brms)
  }))



#Source all files in the auto load funs directory
list.files(file.path(.wd,'analysis/src/funs/auto'),full.names=TRUE) %>%
  walk(source)

conflict_prefer("accumulate", "foreach")
conflict_prefer("ar", "brms")
conflict_prefer("lag", "stats")
conflict_prefer("when", "foreach")


#---- Initialize database ----#
message("Initializing database connection...")

db <- dbConnect(RSQLite::SQLite(), .dbPF, `synchronous` = NULL)

indtb <- tbl(db,'individual') %>%  # Load a tibble with all individual animals
  filter(taxon_canonical_name == "Cervus elaphus") %>% 
  collect()
ind <- indtb %>% 
  pull(individual_id)
# Load the entire event table:

evt0 <- tbl(db, "event_clean") %>%
  filter(individual_id %in% ind) %>% 
  collect()

message("Disconnecting from databse...")
dbDisconnect(db)

ind <- ind[ind  %in% unique(evt0$individual_id)]

yearvec <- c("2019", "2020")

# Add empty columns study:

# log <- read_csv(glue("{.outPF}/niche_log.csv")) #####
ind_det_unscaled_nolog <- list()
ind_det_scaled_nolog <- list()
ind_det_unscaled_log <- list()
ind_det_scaled_log <- list()

for(j in 1:length(unique(ind))){
  
  yr_det_unscaled_nolog <- list()
  yr_det_scaled_nolog <- list()
  yr_det_unscaled_log <- list()
  yr_det_scaled_log <- list()
  
  for(i in 1:2){
    
    #---- Perform analysis ----#
    message(glue("Gathering movement data for individual {ind[j]}, year {i}..."))
    
    scientificname <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(taxon_canonical_name)
    
    studyid <- indtb %>% 
      filter(individual_id == !!ind[j]) %>% 
      pull(study_id)
    
    
    message("Filtering and scaling data...")
    
    # print(i)
    
    # i = 2020
    evt_mod <- evt0 %>% 
      filter(individual_id == ind[j]) %>%
      dplyr::filter(yr == yearvec[i]) %>%
      mutate(timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %T"),
             week = week(timestamp),
             n_indiv_week_year = paste0(individual_id, '_' , week, '_' , yr),
             tmax_scale = scale(tmax),
             tmin_scale = scale(tmin),
             ndvi_scale = scale(ndvi),
             elev_scale = scale(elev)) %>% 
      arrange(timestamp)
    
    
    #-- Fit Multivariate niches ####
    message("Estimating Multivariate niches")
    
    wk <- unique(evt_mod$week)
    # wk = wk[complete.cases(wk)]
    if(length(wk)==0){print(paste0('No data for year ', i, 
                                   ' writing in logfile'))
      
      # tmp_dummy_fail = data.frame(
      #   studyid = studyid, 
      #   individual  = ind[j],
      #   scientificname = scientificname,
      #   year = i,
      #   week = NA,
      #   status = 0)
      # 
      # # logfile_template$studyid <- studyid
      # write.table(tmp_dummy_fail, glue("./out/niche_log.csv"), append = T, 
      #             row.names = F, col.names = F, sep = ",")
      # 
      
    } else { # if no weeks in data
      wk_det_unscaled_nolog <- list()
      wk_det_scaled_nolog <- list()
      wk_det_unscaled_log <- list()
      wk_det_scaled_log <- list()
      
      
      # i <- 10
      for(w in 1:length(wk)){
        # w <- 26
        
        evt_tmp <- evt_mod %>% 
          filter(week == wk[w]) %>%
          select(tmax_scale, tmax, 
                 tmin_scale, tmin,
                 ndvi_scale, ndvi,
                 elev_scale, elev,
                 week, event_id,
                 individual_id, 
                 n_indiv_week_year) %>%
          drop_na(tmax_scale, tmax, tmin_scale, tmin, ndvi_scale, ndvi, 
                  elev_scale, elev)
        
        
        tryCatch({
          
          if(nrow(evt_tmp) > 0){      
            
            wk_det_unscaled_nolog[[w]] <- evt_tmp %>% 
              select(tmax, tmin, ndvi, elev) %>% 
              MVNH_det(log = F) %>% 
              t() %>% 
              as.data.frame() %>% 
              mutate(week = wk[w],
                     individual = unique(evt_tmp$individual_id),
                     year = yearvec[i],
                     scaled = "N",
                     logged = "N")
            
            
            wk_det_scaled_nolog[[w]] <- evt_tmp %>% 
              select(tmax_scale, tmin_scale, ndvi_scale, elev_scale) %>% 
              MVNH_det(log = F) %>% 
              t() %>% 
              as.data.frame() %>% 
              mutate(week = wk[w],
                     individual = unique(evt_tmp$individual_id),
                     year = yearvec[i],
                     scaled = "Y",
                     logged = "N")
            
            wk_det_unscaled_log[[w]] <- evt_tmp %>% 
              select(tmax, tmin, ndvi, elev) %>% 
              MVNH_det(log = T) %>% 
              t() %>% 
              as.data.frame() %>% 
              mutate(week = wk[w],
                     individual = unique(evt_tmp$individual_id),
                     year = yearvec[i],
                     scaled = "N",
                     logged = "T")
            
            wk_det_scaled_log[[w]] <- evt_tmp %>% 
              select(tmax_scale, tmin_scale, ndvi_scale, elev_scale) %>% 
              MVNH_det(log = T) %>% 
              t() %>% 
              as.data.frame() %>% 
              mutate(week = wk[w],
                     individual = unique(evt_tmp$individual_id),
                     year = yearvec[i],
                     scaled = "Y",
                     logged = "Y")
            
            
          } else {# if there's no data
            
            
            print(paste0(unique(evt_tmp$n_indiv_week_year), 
                         ' has NA in niche determinant, moving to the next week'))
          } # else
          
        }, error = function(e){cat(
          glue("ERROR: unspecified error in fitting niche determinant for ind {ind[j]}, yr {i}", 
               "\n"))})
        
      } # for w in wks
      yr_det_unscaled_nolog[[i]] <- do.call("rbind", wk_det_unscaled_nolog)
      yr_det_scaled_nolog[[i]] <- do.call("rbind", wk_det_scaled_nolog)
      yr_det_unscaled_log[[i]] <- do.call("rbind", wk_det_unscaled_log)
      yr_det_scaled_log[[i]] <- do.call("rbind", wk_det_scaled_log)
      
    } # else (if wks > 0)
  } #i (end loop through years)
  ind_det_unscaled_nolog[[j]] <- do.call("rbind", yr_det_unscaled_nolog)
  ind_det_scaled_nolog[[j]] <- do.call("rbind", yr_det_scaled_nolog)
  ind_det_unscaled_log[[j]] <- do.call("rbind", yr_det_unscaled_log)
  ind_det_scaled_log[[j]] <- do.call("rbind", yr_det_scaled_log)
  
}#j (end loop through individuals)

det_unscaled_nolog <- do.call("rbind", ind_det_unscaled_nolog)
det_scaled_nolog <- do.call("rbind", ind_det_scaled_nolog)
det_unscaled_log <- do.call("rbind", ind_det_unscaled_log)
det_scaled_log <- do.call("rbind", ind_det_scaled_log)



dbbmms <- read_csv(.varPF) %>% 
  mutate(species = case_when( # correct species names
    study_id == 1442516400 ~ "Anser caerulescens",
    study_id == 1233029719 ~ "Odocoileus virginianus",
    study_id == 1631574074 ~ "Ursus americanus",
    study_id == 1418296656 ~ "Numenius americanus",
    study_id == 474651680  ~ "Odocoileus virginianus",
    study_id == 1044238185 ~ "Alces alces",
    TRUE ~ species
  ))%>% 
  mutate(species = case_when(
    species == "Chen caerulescens" ~ "Anser caerulescens",
    TRUE ~ species
  )) %>% 
  mutate(year_f = as.character(year))


breadth_unscaled_nolog <- det_unscaled_nolog %>%
  mutate(tmax_mvnh = tmax,
         tmin_mvnh = tmin,
         elev_mvnh = elev,
         ndvi_mvnh = ndvi) %>%
  select(!c(tmax, tmin, elev, ndvi)) %>%
  filter(studyid != 351564596) %>%
  filter(studyid != 1891587670) %>%
  mutate(ind_f = as.factor(individual)) %>%  # create factor version of ind for REs)
  left_join(dbbmms, by = c("individual" = "ind_id", 
                           "year" = "year_f", 
                           "week" = "wk")) %>% 
  filter(!is.infinite(total)) %>% 
  mutate(
    # sqrt_breadth = sqrt(total), #get sqrt weekly niche breadth
    breadth = total,
    breadth_scale = scale(total), # standardize it
    sg_norm = scale(sg / cbg_area), # normalize safegraph data by size of the CBG
    # log_sg_norm = log(sg_norm),
    ghm_scale = scale(ghm),
    ndvi_scale = scale(ndvi),
    lst_scale = scale(lst),
    tmax_scale = scale(tmax),
    tmin_scale = scale(tmin),
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


breadth_scaled_nolog <- det_scaled_nolog %>%
  mutate(tmax_mvnh = tmax_scale,
         tmin_mvnh = tmin_scale,
         elev_mvnh = elev_scale,
         ndvi_mvnh = ndvi_scale) %>%
  select(!c(tmax_scale, tmin_scale, elev_scale, ndvi_scale)) %>%
  filter(studyid != 351564596) %>%
  filter(studyid != 1891587670) %>%
  mutate(ind_f = as.factor(individual)) %>%  # create factor version of ind for REs)
  left_join(dbbmms, by = c("individual" = "ind_id", 
                           "year" = "year_f", 
                           "week" = "wk")) %>% 
  filter(!is.infinite(total)) %>% 
  mutate(
    # sqrt_breadth = sqrt(total), #get sqrt weekly niche breadth
    breadth = total,
    breadth_scale = scale(total), # standardize it
    sg_norm = scale(sg / cbg_area), # normalize safegraph data by size of the CBG
    # log_sg_norm = log(sg_norm),
    ghm_scale = scale(ghm),
    ndvi_scale = scale(ndvi),
    lst_scale = scale(lst),
    tmax_scale = scale(tmax),
    tmin_scale = scale(tmin),
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


breadth_unscaled_log <- det_unscaled_log %>%
  mutate(tmax_mvnh = tmax,
         tmin_mvnh = tmin,
         elev_mvnh = elev,
         ndvi_mvnh = ndvi) %>%
  select(!c(tmax, tmin, elev, ndvi)) %>%
  filter(studyid != 351564596) %>%
  filter(studyid != 1891587670) %>%
  mutate(ind_f = as.factor(individual)) %>%  # create factor version of ind for REs)
  left_join(dbbmms, by = c("individual" = "ind_id", 
                           "year" = "year_f", 
                           "week" = "wk")) %>% 
  filter(!is.infinite(total)) %>% 
  mutate(
    # sqrt_breadth = sqrt(total), #get sqrt weekly niche breadth
    breadth = total,
    breadth_scale = scale(total), # standardize it
    sg_norm = scale(sg / cbg_area), # normalize safegraph data by size of the CBG
    # log_sg_norm = log(sg_norm),
    ghm_scale = scale(ghm),
    ndvi_scale = scale(ndvi),
    lst_scale = scale(lst),
    tmax_scale = scale(tmax),
    tmin_scale = scale(tmin),
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


breadth_scaled_log <- det_scaled_log %>%
  mutate(tmax_mvnh = tmax_scale,
         tmin_mvnh = tmin_scale,
         elev_mvnh = elev_scale,
         ndvi_mvnh = ndvi_scale) %>%
  select(!c(tmax_scale, tmin_scale, elev_scale, ndvi_scale)) %>%
  filter(studyid != 351564596) %>%
  filter(studyid != 1891587670) %>%
  mutate(ind_f = as.factor(individual)) %>%  # create factor version of ind for REs)
  left_join(dbbmms, by = c("individual" = "ind_id", 
                           "year" = "year_f", 
                           "week" = "wk")) %>% 
  filter(!is.infinite(total)) %>% 
  mutate(
    # sqrt_breadth = sqrt(total), #get sqrt weekly niche breadth
    breadth = total,
    breadth_scale = scale(total), # standardize it
    sg_norm = scale(sg / cbg_area), # normalize safegraph data by size of the CBG
    # log_sg_norm = log(sg_norm),
    ghm_scale = scale(ghm),
    ndvi_scale = scale(ndvi),
    lst_scale = scale(lst),
    tmax_scale = scale(tmax),
    tmin_scale = scale(tmin),
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


#---- Models ----#
conflict_prefer("ar", "brms")
conflict_prefer("collapse", "dplyr")
conflict_prefer("lag", "stats")
conflict_prefer("Position", "base")

form_s <-  bf(breadth_scale ~ sg_norm + ghm_scale + (1 |grp) + ar(time = week, gr = grp))
# form_n <-  bf(breadth ~ sg_norm + ghm_scale + (1 |grp) + ar(time = week, gr = grp))

  

# fit model
mod_us_nl_s <- brm(
  form_s,
  data = breadth_unscaled_nolog,
  # family = Gamma(link = "log"),
  inits = 0,
  cores = 1,
  iter = .iter,
  thin = .thin
)

mod_s_nl_s <- brm(
  form_s,
  data = breadth_scaled_nolog,
  # family = Gamma(link = "log"),
  inits = 0,
  cores = 1,
  iter = .iter,
  thin = .thin
)

mod_us_l_s <- brm(
  form_s,
  data = breadth_unscaled_log,
  # family = Gamma(link = "log"),
  inits = 0,
  cores = 1,
  iter = .iter,
  thin = .thin
)

mod_s_l_s <- brm(
  form_s,
  data = breadth_scaled_log,
  # family = Gamma(link = "log"),
  inits = 0,
  cores = 2,
  iter = .iter,
  thin = .thin
)

fe_us_nl_s <- fixef(mod_us_nl_s) #get fixed effects
df_us_nl_s <- tibble(# SG OUT
                "add_sg_norm"=as.numeric(fe_us_nl_s["sg_norm", "Estimate"]),
                "add_sg_norm_lci"=fe_us_nl_s["sg_norm", "Q2.5"],
                "add_sg_norm_uci"=fe_us_nl_s["sg_norm", "Q97.5"],
                
                # GHM OUT
                "add_ghm_scale"=as.numeric(fe_us_nl_s["ghm_scale", "Estimate"]),
                "add_ghm_scale_lci"=fe_us_nl_s["ghm_scale", "Q2.5"],
                "add_ghm_scale_uci"=fe_us_nl_s["ghm_scale", "Q97.5"]) %>% 
  mutate(add_sg_sig = case_when(
      (add_sg_norm_lci < 0 & 0 < add_sg_norm_uci) ~ "N",
      TRUE ~ "Y"),
    add_ghm_sig = case_when(
      (add_ghm_scale_lci < 0 & 0 < add_ghm_scale_uci) ~ "N",
      TRUE ~ "Y"),
    model = "us_nl_s")

fe_s_nl_s <- fixef(mod_s_nl_s) #get fixed effects
df_s_nl_s <- tibble(# SG OUT
  "add_sg_norm"=as.numeric(fe_s_nl_s["sg_norm", "Estimate"]),
  "add_sg_norm_lci"=fe_s_nl_s["sg_norm", "Q2.5"],
  "add_sg_norm_uci"=fe_s_nl_s["sg_norm", "Q97.5"],
  
  # GHM OUT
  "add_ghm_scale"=as.numeric(fe_s_nl_s["ghm_scale", "Estimate"]),
  "add_ghm_scale_lci"=fe_s_nl_s["ghm_scale", "Q2.5"],
  "add_ghm_scale_uci"=fe_s_nl_s["ghm_scale", "Q97.5"]) %>% 
  mutate(add_sg_sig = case_when(
    (add_sg_norm_lci < 0 & 0 < add_sg_norm_uci) ~ "N",
    TRUE ~ "Y"),
    add_ghm_sig = case_when(
      (add_ghm_scale_lci < 0 & 0 < add_ghm_scale_uci) ~ "N",
      TRUE ~ "Y"),
    model = "s_nl_s")

fe_us_l_s <- fixef(mod_us_l_s) #get fixed effects
df_us_l_s <- tibble(# SG OUT
  "add_sg_norm"=as.numeric(fe_us_l_s["sg_norm", "Estimate"]),
  "add_sg_norm_lci"=fe_us_l_s["sg_norm", "Q2.5"],
  "add_sg_norm_uci"=fe_us_l_s["sg_norm", "Q97.5"],
  
  # GHM OUT
  "add_ghm_scale"=as.numeric(fe_us_l_s["ghm_scale", "Estimate"]),
  "add_ghm_scale_lci"=fe_us_l_s["ghm_scale", "Q2.5"],
  "add_ghm_scale_uci"=fe_us_l_s["ghm_scale", "Q97.5"]) %>% 
  mutate(add_sg_sig = case_when(
    (add_sg_norm_lci < 0 & 0 < add_sg_norm_uci) ~ "N",
    TRUE ~ "Y"),
    add_ghm_sig = case_when(
      (add_ghm_scale_lci < 0 & 0 < add_ghm_scale_uci) ~ "N",
      TRUE ~ "Y"),
    model = "us_l_s")

fe_s_l_s <- fixef(mod_s_l_s) #get fixed effects
df_s_l_s <- tibble(# SG OUT
  "add_sg_norm"=as.numeric(fe_s_l_s["sg_norm", "Estimate"]),
  "add_sg_norm_lci"=fe_s_l_s["sg_norm", "Q2.5"],
  "add_sg_norm_uci"=fe_s_l_s["sg_norm", "Q97.5"],
  
  # GHM OUT
  "add_ghm_scale"=as.numeric(fe_s_l_s["ghm_scale", "Estimate"]),
  "add_ghm_scale_lci"=fe_s_l_s["ghm_scale", "Q2.5"],
  "add_ghm_scale_uci"=fe_s_l_s["ghm_scale", "Q97.5"]) %>% 
  mutate(add_sg_sig = case_when(
    (add_sg_norm_lci < 0 & 0 < add_sg_norm_uci) ~ "N",
    TRUE ~ "Y"),
    add_ghm_sig = case_when(
      (add_ghm_scale_lci < 0 & 0 < add_ghm_scale_uci) ~ "N",
      TRUE ~ "Y"),
    model = "s_l_s")

fe_run <- fixef(out$model) #get fixed effects
df_run <- tibble(# SG OUT
  "add_sg_norm"=as.numeric(fe_run["sg_norm", "Estimate"]),
  "add_sg_norm_lci"=fe_run["sg_norm", "Q2.5"],
  "add_sg_norm_uci"=fe_run["sg_norm", "Q97.5"],
  
  # GHM OUT
  "add_ghm_scale"=as.numeric(fe_run["ghm_scale", "Estimate"]),
  "add_ghm_scale_lci"=fe_run["ghm_scale", "Q2.5"],
  "add_ghm_scale_uci"=fe_run["ghm_scale", "Q97.5"]) %>% 
  mutate(add_sg_sig = case_when(
    (add_sg_norm_lci < 0 & 0 < add_sg_norm_uci) ~ "N",
    TRUE ~ "Y"),
    add_ghm_sig = case_when(
      (add_ghm_scale_lci < 0 & 0 < add_ghm_scale_uci) ~ "N",
      TRUE ~ "Y"),
    model = "run")

(comb <- rbind(df_run, df_us_l_s, df_us_nl_s, df_s_nl_s, df_s_l_s))

#-- Check Fits --#
pp_check(out$model)
pp_check(mod_s_l_s)

#---- Finalize script ----#

message(glue('Script complete in {diffmin(t0)} minutes'))
# load("out/fixniche.rdata")
