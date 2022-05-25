#!/bin/bash

#SBATCH -t 8:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition pi_jetz
#SBATCH -c 10
#SBATCH --mem-per-cpu 10G
#SBATCH -J _sigle_species-Loop_2022-05-16

# Load conda env
module load miniconda
conda activate brms

# Declare WD
wd=~/project/covid-19_movement

#ensure we're in wd
cd $wd

# Execute calc size script/
Rscript analysis/src/single_species_loop_niche.r
