#!/bin/bash

#SBATCH -t 01:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH -c 1
#SBATCH --mem-per-cpu 50G
#SBATCH -J prep_covid_2022-01-20

# Loda conda env
module load miniconda
conda activate spatial

# Execute cleaning script (non parallel)
Rscript $wd/analysis/src/01-prep_and_clean.r --db $wd/processed_data/mosey_mod.db