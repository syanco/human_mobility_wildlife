#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition day
#SBATCH -c 24
#SBATCH --mem-per-cpu 10G
#SBATCH -J fit_space_use_mods_INT

# Load conda env
module load miniconda
conda activate brms

# Declare WD
wd=~/project/covid-19_movement

cd $wd

# Execute model script 
Rscript $wd/analysis/src/workflow/fit-space-use-interactive-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area_interactive 24 5 10000 5
