require(tidyverse)
require(stringi)

setwd('/gpfs/loomis/pi/jetz/de293/Anthropause/out/NVMH_subsampling/')

all_niches_subsampled <- read.csv('all_niches_subsampled.csv')

names(all_niches_subsampled) <- c('RowId', 'total', 'tmax_scale', 'tmin_scale', 'ndvi_scale',
                                  'elev_scale', 'cor', 'week', 'individual', 'year', 'iteration',
                                  'n_sample', 'n_indiv_week_year')

require(plyr)

# ddply median niche breadth per iteration?
# save plot to file without using ggsave
p <-
  all_niches_subsampled %>% 
  ggplot(aes(y = (total), x = n_sample)) +
  geom_boxplot(aes(fill = factor(n_sample)), outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +theme_bw() +
  xlab('Sample size')+
  ylab(('Niche breadth'))+
  ggtitle('Niche breadth across sample size') + 
  theme(plot.title = element_text(size = 18, face = "italic"),
        element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(face = 'bold'),
        axis.title.x = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold')) +
  theme(legend.position = 'none') + xlim(10, 50)
# xlim(10, 30, 50)
# + scale_y_continuous(breaks = seq(0, 100, by = 20))
png('/gpfs/loomis/pi/jetz/de293/Anthropause/out/Figures/All_niches_2.png')
print(p)
dev.off()


all_niches_subsampled %>% 
  ggplot(aes(y = (total), x = n_sample)) +
  geom_boxplot(aes(fill = factor(n_sample)), outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE) +theme_bw() +
  xlab('Sample size')+
  ylab(('Niche breadth'))+
  ggtitle('Niche breadth across sample size') + 
  theme(plot.title = element_text(size = 18, face = "italic"),
        element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(face = 'bold'),
        axis.title.x = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold')) +
  theme(legend.position = 'none') + xlim(10, 50)


ggsave('/gpfs/loomis/pi/jetz/de293/Anthropause/out/Figures/All_niches.png'
       , width = 6
       , height = 7
       , dpi = 600
)
