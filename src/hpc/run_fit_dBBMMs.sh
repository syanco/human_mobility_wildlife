#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition pi_jetz
#SBATCH -c 24
#SBATCH --mem-per-cpu 10G
#SBATCH -J fit_dBBMMs

# Load conda env
module load miniconda
conda activate covid

# Declare WD
wd=~/project/covid-19_movement

#copy db to tmp
cp $wd/processed_data/mosey_mod.db /tmp/

# Execute cleaning script (non parallel)
Rscript $wd/analysis/src/workflow/fit-dBBMMs.r /tmp/mosey_mod.db $wd/out 24
