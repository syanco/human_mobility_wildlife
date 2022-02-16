#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition pi_jetz
#SBATCH -c 4
#SBATCH --mem-per-cpu 10G
#SBATCH -J quick_mod_2022-02-16

# Load conda env
module load miniconda
conda activate brms

# Declare WD
wd=~/project/covid-19_movement

#ensure we're in wd
cd $wd

# Execute calc size script/
Rscript analysis/src/quick_model.r

