#!/bin/bash

#SBATCH --job-name=filter_data_mins
#SBATCH -t 05:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu
#SBATCH -c 1
##SBATCH --mem-per-cpu 50G            # add this back in if needed

module load R/4.1.3 gdal/2.2.3 proj/5.2

# Declare WD
export wd=/scratch/julietcohen/covid_movement/human_mobility_wildlife
export src=$wd/src/workflow

cd $wd

# Execute cleaning script (non parallel)
Rscript $src/filter_data_mins.R --db /tmp/mosey_mod.db 30 5
