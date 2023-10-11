#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition day
#SBATCH -c 6
#SBATCH --mem-per-cpu 10G
#SBATCH -J fit_space_use_mods_DOT

# Load conda env
module load miniconda
conda activate brms

# Declare WD
wd=~/project/covid-19_movement

cd $wd

# Execute model script 
Rscript $wd/analysis/src/workflow/refit-space-use-dot-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area_dot 24 5 10000 5
