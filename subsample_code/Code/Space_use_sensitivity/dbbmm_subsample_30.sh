#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user diego.ellissoto@yale.edu
#SBATCH --partition day
#SBATCH -c 15
#SBATCH --mem-per-cpu 20G
#SBATCH -J dbbmm_sensitivity_week_30_sub

# Load conda env
module load miniconda
conda activate covid

Rscript dbbmm_subsample_30.R
