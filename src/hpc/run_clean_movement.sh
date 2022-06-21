#!/bin/bash

#SBATCH -t 01:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH -c 1
#SBATCH --mem-per-cpu 50G

# Load conda env
module load miniconda
conda activate spatial

# Declare WD
wd=~/project/covid-19_movement

# Execute cleaning script (non parallel)
Rscript $wd/analysis/src/workflow/clean_movement.r --db $wd/processed_data/mosey_mod.db