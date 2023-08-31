#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user diego.ellissoto@yale.edu
#SBATCH --partition day
#SBATCH -c 1
#SBATCH --mem-per-cpu 20G
SBATCH -J write_indiv_week_year

# Load conda env
module load miniconda
conda activate covid

Rscript write_indiv_week_year.R
