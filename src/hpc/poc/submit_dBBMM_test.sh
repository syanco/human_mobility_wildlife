#!/bin/bash

#SBATCH -t 10:00:00
#SBATCH --mail-type None
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition bigmem
#SBATCH -c 24
#SBATCH --mem-per-cpu=20G
#SBATCH -J fit_dBBMMs_tjl

# Load conda env
module load miniconda
conda activate covid

# Declare WD
wd=~/project/covid-19_movement

# Execute cleaning script (non parallel)
Rscript $wd/analysis/src/0X-fit_dBBMMs.r $wd/processed_data/mosey_mod.db $wd/out 24
