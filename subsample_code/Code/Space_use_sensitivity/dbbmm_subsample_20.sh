#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user diego.ellissoto@yale.edu
#SBATCH --partition day
#SBATCH -c 15
#SBATCH --mem-per-cpu 20G
#SBATCH -J dbbmm_sensitivity

# Load conda env
module load miniconda
conda activate covid

Rscript dbbmm_subsample_20.R
