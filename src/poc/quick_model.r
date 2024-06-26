library(brms)
library(tidyverse)
library(ape)
library(lubridate)

vertnet_tax <- read_csv("raw_data/vertlife_taxonomies.csv")

size <- read_csv("out/dbbmm_size.csv") %>% 
  mutate(
    # trt_new = gsub('_.*','',trt),
    year_f = factor(year),
    # trt_new = fct_relevel(trt_new, "pre.ld", "ld", "post.ld")
  )

# join species in dataset to the taxonomy
foc_sp <- vertnet_tax %>% 
  inner_join(size, by = c("scientificname" = "species"))

#read in trees
birds <- read.nexus("raw_data/tree-pruner-05469463-d432-42ac-a726-58ff39a822b8/output.nex")
birds <- birds$tree_4345
mammals <- read.nexus("raw_data/tree-pruner-baceaa54-30cf-4bdd-a253-a20897aa6496/output.nex")
mammals <- mammals$tree_2272
combphylo <- bind.tree(birds, mammals)
phylo_vcov <- ape::vcv.phylo(combphylo)

# fix species names to match phylogeny
size <- size %>% 
  mutate(sp2 = gsub(" ", "_", species)) %>% 
  # also remove any species not in tree (i.e. NAs)
  # filter(sp2 %in% rownames(phylo_vcov)) %>% 
  mutate(ind_f = as.factor(ind_id),
         wk_n = as.numeric(substring(wk,2)),
         ts = parse_date_time(paste(year, wk, 01, sep = "-"), "%Y-%U-%u" )) %>%
  distinct()


#- Load ain the species trait data
traits <- read_csv("raw_data/anthropause_data_sheet.csv") %>% 
  mutate(mig_mod = case_when(migratory == "Partial" ~ "migratory",
                             .$migratory == "non-migratory" ~ "non-migratory",
                             .$migratory == "Complete" ~ "migratory",
                             .$migratory == "migratory" ~ "migratory",
                             .$migratory == "Semi-nomadic" ~ "non-migratory",
                             .$migratory == "Unknown" ~ NA_character_))

size <- size %>% 
  left_join(traits, by = c("species" = "Species"))

size_res <- size %>%
  filter(mig_mod == "non-migratory") %>%
  filter(study_id != 351564596) %>%
  filter(study_id != 1891587670) %>%
  mutate(log_area = log(area),
         sg_norm = area/cbg_area,
         log_sg_norm = log(sg_norm),
         grp = as.factor(paste(ind_f, year, sep = "_")),
         year_f = as.factor(year)) %>%
  group_by(grp) %>%
  arrange(.by_group = T)

# form <- bf(area ~ year_f + s(wk_n, by = year_f) + (1|ind_f) + (1|gr(sp2, cov = phylo_vcov)))
# form <- bf(area ~ pop*sg + ndvi + year_f + (1|sp2) + (1|ind_f) + ar(time = ts, p = 1, gr = ind_f:year_f, cov = F))
form <- bf(log_area ~ sg_norm*sp2 + ghm + year_f + ndvi + lst + (1|grp) + ar(time = wk, gr=grp)) 
# form <- bf(log_sg_norm ~ 1 + (1|ind_f) + ar(time = ts ,gr=ind_f)) 

           # + ar(time = ts, gr = ind_f, cov = F))

mod_res <- brm(form, 
           # data2 = list(phylo_vcov = phylo_vcov),
           data = size_res,
           # family = Gamma(link = "log"),
           inits = 0,
           cores = 4,
           iter = 500)

pp_check(mod_res)
plot(mod_res)

save("out/quick_mod.rdata")
