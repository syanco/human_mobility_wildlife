#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user diego.ellissoto@yale.edu
#SBATCH --partition day
#SBATCH -c 1
#SBATCH --mem-per-cpu 10G
SBATCH -J one_p_p_day

# Load conda env
module load miniconda
conda activate covid

Rscript 2023_one_p_p_day.R
