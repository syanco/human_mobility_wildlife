require(tidyverse)
require(stringi)
require(plyr)

setwd('/gpfs/loomis/pi/jetz/de293/Anthropause/out/NVMH_subsampling/')

all_niches_subsampled <- read.csv('all_niches_subsampled.csv')

names(all_niches_subsampled) <- c('RowId', 'total', 'tmax_scale', 'tmin_scale', 'ndvi_scale',
                                  'elev_scale', 'cor', 'week', 'individual', 'year', 'iteration',
                                  'n_sample', 'n_indiv_week_year')

write.csv( sample_n(all_niches_subsampled, 100000), file = '/gpfs/loomis/pi/jetz/de293/Anthropause/out/NVMH_subsampling/all_niches_subsampled_sub_100000.csv')
           
