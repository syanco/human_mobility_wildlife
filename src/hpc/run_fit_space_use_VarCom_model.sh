#!/bin/bash

#SBATCH -t 2-
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition pi_jetz,day
#SBATCH -c 4
#SBATCH --mem-per-cpu 10G
#SBATCH -J fit_space_use_trait_mod

# Load conda env
module load miniconda
conda activate brms

# Declare WD
wd=/gpfs/loomis/pi/jetz/sy522/covid-19_movement

cd $wd

# Execute model script 
Rscript $wd/analysis/src/workflow/fit-space-use-VarCom-model.r $wd/out/dbbmm_size.csv $wd/out/single_species_models $wd/raw_data/anthropause_data_sheet.csv 4 4000 2
