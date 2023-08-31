#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user diego.ellissoto@yale.edu
#SBATCH --partition day
#SBATCH -c 20
#SBATCH --mem-per-cpu 10G
#SBATCH -J niche_subsample_20

# Load conda env
module load miniconda
conda activate covid

Rscript niche_subsample_20.R
