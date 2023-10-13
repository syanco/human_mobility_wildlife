#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition day
#SBATCH -c 8
#SBATCH --mem-per-cpu 10G
#SBATCH -J refit_add_niche

# Load conda env
module load miniconda
conda activate brms

# Declare WD
wd=~/project/covid-19_movement

cd $wd

# Execute model script 
Rscript $wd/analysis/src/workflow/refit-niche-breadth-additive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_additive 8 
