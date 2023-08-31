#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user diego.ellissoto@yale.edu
#SBATCH --partition day
#SBATCH -c 20
#SBATCH --mem-per-cpu 10G
#SBATCH -J anthro_subs_points_2023_50p

# Load conda env
module load miniconda
conda activate covid

Rscript anthro_point_subsampling.R
