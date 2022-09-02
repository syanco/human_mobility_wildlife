#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition pi_jetz,day
#SBATCH -c 24
#SBATCH --mem-per-cpu 10G
#SBATCH -J fit_niche_bredth_add_mods

# Load conda env
module load miniconda
conda activate brms

# Declare WD
wd=/gpfs/loomis/pi/jetz/sy522/covid-19_movement

cd $wd

# Execute model script 
Rscript $wd/analysis/src/workflow/fit-niche-breadth-additive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_controlled 24 10 10000 5
