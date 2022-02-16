library(brms)
library(tidyverse)
library(ape)

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

#print species by phylum
foc_sp %>% 
  filter(group == "Birds") %>% 
  select(scientificname) %>% 
  unique() %>%
  print(n=25)

foc_sp %>% 
  filter(group == "Mammals") %>% 
  select(scientificname) %>% 
  unique() %>%
  print(n=25)

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
  filter(sp2 %in% rownames(phylo_vcov)) %>% 
  mutate(ind_f = as.factor(ind_id),
         wk_n = as.numeric(substring(wk,2)))

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

size2 <- size %>% 
  filter(mig_mod == "migratory") 
# form1 <- bf(area ~ trt_new*year_f + phen + (1|gr(sp2, cov = phylo_vcov)),
#             # ndt ~ trt_new*year_f + phen + (1|gr(sp2, cov = phylo_vcov)))
#             ndt ~ 1)
# 
# #use this to get stan code (to extract param names for inits)
# make_stancode(form1, 
#               data = size, 
#               data2 = list(phylo_vcov = phylo_vcov), 
#               family = shifted_lognormal())
# 
# inits <- list(sigma = 1 
#               # sd_1 = 1,
#               # sd_2 = 1
#               )
# 
# inits_list <- list(inits, inits, inits, inits)
# 
# mod1 <- brm(
#   # area ~ trt_new*year_f + (1|gr(sp2, cov = phylo_vcov)), 
#   form1,
#   data = size, 
#   # inits = inits_list,
#   family = shifted_lognormal(),
#   # family = gaussian(),
#   data2 = list(phylo_vcov = phylo_vcov),
#   cores = 2,
#   control = list(adapt_delta = 0.95))
# 
# mod1
# conditional_effects(mod1, points = T)
# pp_check(mod1, type = "stat", stat = "mean")
# pp_check(mod1)

form2 <- bf(area ~ year_f + s(wk_n, by = year_f) + (1|gr(sp2, cov = phylo_vcov)))

mod2 <- brm(form2, 
            data2 = list(phylo_vcov = phylo_vcov),
            data = size,
            family = Gamma(link = "log"),
            inits = 0,
            cores = 4)

mod2
conditional_effects(mod2, points = T)
pp_check(mod2, type = "stat", stat = "mean")
pp_check(mod2)

# let slpes vary as random effect of phylogeny
form3 <- bf(area ~ trt_new*year_f*mig_mod + phen + (1+trt_new*year_f*mig_mod|gr(sp2, cov = phylo_vcov))
            + (1|ind_f))

mod3 <- brm(form3, 
            data2 = list(phylo_vcov = phylo_vcov),
            data = size,
            family = Gamma(link = "log"),
            inits = 0,
            cores = 2)

mod3
conditional_effects(mod3, points = T)
pp_check(mod3, type = "stat", stat = "mean")
pp_check(mod3)



###############################
x <- size %>% 
  group_by(species) %>% 
  summarize(n = n_distinct(ind_f))
