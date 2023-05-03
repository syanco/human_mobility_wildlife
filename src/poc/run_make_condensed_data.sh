#!/bin/bash

#SBATCH -t 2:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH -c 1
#SBATCH --mem-per-cpu 25G

# Load conda env
module load miniconda
conda activate covid

# Declare WD
wd=~/project/covid-19_movement
src=$wd/analysis/src/poc

# Move to WD
cd $wd

# Execute condense script (non parallel)
Rscript $src/make_condensed_data.r 