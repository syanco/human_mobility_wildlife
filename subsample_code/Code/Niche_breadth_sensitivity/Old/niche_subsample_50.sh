#!/bin/bash

#SBATCH -t 2-
#SBATCH --mail-type ALL
#SBATCH --mail-user diego.ellissoto@yale.edu
#SBATCH --partition pi_jetz
#SBATCH -c 20
#SBATCH --mem-per-cpu 5G
#SBATCH -J nvmh_sensitivity

# Load conda env
module load miniconda
conda activate covid

Rscript niche_subsample_50.R
