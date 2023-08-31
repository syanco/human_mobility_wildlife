#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user diego.ellissoto@yale.edu
#SBATCH --partition pi_jetz
#SBATCH -c 1
#SBATCH --mem-per-cpu 120G
#SBATCH -J boxplot_niche_breadth

# Load conda env
module load miniconda
conda activate covid

Rscript boxplot_niche_breadth_subsampling.R
