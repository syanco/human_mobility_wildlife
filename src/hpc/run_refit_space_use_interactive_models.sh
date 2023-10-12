#!/bin/bash

#SBATCH -t 12:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition day
#SBATCH -c 8
#SBATCH --mem-per-cpu 10G
#SBATCH -J fit_space_use_mods_ADD

# Load conda env
module load miniconda
conda activate brms

# Declare WD
wd=~/project/covid-19_movement

cd $wd

# Execute model script 
Rscript $wd/analysis/src/workflow/refit-space-use-interactive-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area_additive 8
