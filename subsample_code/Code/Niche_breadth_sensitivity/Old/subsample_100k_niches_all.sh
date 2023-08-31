#!/bin/bash

#SBATCH -t 06:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user diego.ellissoto@yale.edu
#SBATCH --partition pi_jetz
#SBATCH -c 1
#SBATCH --mem-per-cpu 50G
#SBATCH -J subsample_niches

# Load conda env
module load miniconda
conda activate covid

Rscript subsample_100k_niches_all.R
