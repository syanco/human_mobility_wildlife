#!/bin/bash

#SBATCH -t 3-
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition pi_jetz
#SBATCH -c 24
#SBATCH --mem-per-cpu 10G
#SBATCH -J fit_dBBMMs_2022-01-28

# Load conda env
module load miniconda
conda activate covid

# Declare WD
wd=~/project/covid-19_movement

# Execute cleaning script (non parallel)
Rscript $wd/analysis/src/0X-fit_dBBMMs.r $wd/processed_data/mosey_mod.db $wd/out 24
