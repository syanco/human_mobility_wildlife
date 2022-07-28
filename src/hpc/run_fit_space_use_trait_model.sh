#!/bin/bash

#SBATCH -t 6:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition day
#SBATCH -c 4
#SBATCH --mem-per-cpu 10G
#SBATCH -J fit_space_use_trait_mod

# Load conda env
module load miniconda
conda activate brms

# Declare WD
wd=~/project/covid-19_movement

cd $wd

# Execute model script 
Rscript $wd/analysis/src/workflow/fit-space-use-trait-model.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area $wd/raw_data/anthropause_data_sheet.csv 4 10000 5
