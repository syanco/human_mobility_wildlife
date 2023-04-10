#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition day
#SBATCH -c 24
#SBATCH --mem-per-cpu 10G
#SBATCH -J fit_niche_bredth_dot_mods

# Load conda env
module load miniconda
conda activate brms

# Declare WD
wd=/gpfs/loomis/pi/jetz/sy522/covid-19_movement

cd $wd

# Execute model script 
Rscript $wd/analysis/src/workflow/fit-niche-breadth-dot-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_dot 24 5 10000 5
