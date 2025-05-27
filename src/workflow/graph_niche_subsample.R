#---- Graph Niche Subsample ----#

library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)

niche_sub_10 <- read_csv("out/niche_subsamples/niche_determinant_anthropause_subsample_10.csv")
niche_sub_20 <- read_csv("out/niche_subsamples/niche_determinant_anthropause_subsample_20.csv")
niche_sub_30 <- read_csv("out/niche_subsamples/niche_determinant_anthropause_subsample_30.csv")
niche_sub_40 <- read_csv("out/niche_subsamples/niche_determinant_anthropause_subsample_40.csv")
niche_sub_50 <- read_csv("out/niche_subsamples/niche_determinant_anthropause_subsample_50.csv")

# combine all into one df
niche_sub <- rbind(niche_sub_10, niche_sub_20, niche_sub_30, niche_sub_40, niche_sub_50)

# pivot to wide form and calc individual and week specific diffs across subsamples
sub_wide <- niche_sub %>% 
  pivot_wider(id_cols = c(studyid, individual, year, week), 
              names_from = n, 
              names_prefix = "b", 
              values_from = total) %>% 
  mutate(comp_40 = b50-b40,
         comp_30 = b50-b40,
         comp_20 = b50-b20,
         comp_10 = b50-b10)

# pivot back to long for graphing
sub_plot_df <- sub_wide %>% 
  pivot_longer(cols = c(comp_40, comp_30, comp_20, comp_10), names_to = "subsample") %>% 
  mutate(subsample = fct_rev(subsample))

ggplot(sub_plot_df)+
  geom_boxplot(aes(y = log(value+0.00000000001), x = subsample))

m1 <- lm(log(value+0.00000000001) ~ subsample, data = sub_plot_df)


summary(m1)


sub_plot <- sub_wide %>% 
  pivot_longer(cols = c(b50, b40, b30, b20, b10), names_to = "subsample") %>% 
  mutate(subsample = fct_rev(subsample))


ggplot(sub_plot)+
  geom_boxplot(aes(y = log(value+0.00000000001), x = subsample))

size_mod <- lmerTest::lmer(log(value+0.00000000001) ~ 0 + subsample + (1|individual/year/week), data = sub_plot)

(mod_sum <- summary(size_mod))

size_ci <- confint(size_mod)

size_mod_out <- as.data.frame(size_ci)[5:9,] %>% 
  mutate(est = mod_sum$coefficients[,1]) %>% 
  rownames_to_column(var = "sample") %>% 
  mutate(x_lab = case_when(sample == "subsampleb10" ~ "subsample 10",
                           sample == "subsampleb20" ~ "subsample 20",
                           sample == "subsampleb30" ~ "subsample 30",
                           sample == "subsampleb40" ~ "subsample 40",
                           sample == "subsampleb50" ~ "subsample 50"))

ggplot(size_mod_out)+
  geom_point(aes(x = x_lab, y = est), size = 3) +
  geom_errorbar(aes(x = x_lab, ymin = `2.5 %`, ymax = `97.5 %`), width = 0.3, size = 1.5)+
  ylab("log(niche breadth)")+
  xlab("")+
  theme_classic()

ggsave("out/niche_subsample_plot.png", width = 9, height = 6)
