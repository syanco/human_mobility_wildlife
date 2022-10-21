#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition pi_jetz,day
#SBATCH -c 4
#SBATCH --mem-per-cpu 10G
#SBATCH -J fit-intra-ind-int-niche

# Load conda env
module load miniconda
conda activate brms

# Declare WD
wd=/gpfs/loomis/pi/jetz/sy522/covid-19_movement

cd $wd

# Execute model script 
Rscript $wd/analysis/src/workflow/fit_intra_ind_mod_interactive_niche.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 4 10000 5