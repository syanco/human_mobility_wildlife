#!/bin/bash

#SBATCH --job-name=bckgrnd-pts
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scott.yanco@yale.edu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000g -t 2-
#SBATCH --partition=pi_jetz

module load miniconda
conda activate covid

# Declare WD
wd=~/project/covid-19_movement
src=$wd/analysis/src/workflow

# Move to WD
cd $wd

Rscript $src/generate-background-points.r