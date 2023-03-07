#!/bin/bash

#SBATCH -t 2:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition day
#SBATCH -c 1
#SBATCH --mem-per-cpu 100G

# module load R/4.1.0-foss-2020b
module load miniconda
conda activate covid

# Declare WD
wd=~/project/covid-19_movement
src=$wd/analysis/src/workflow

# Move to WD
cd $wd

Rscript $src/compute-cbg-area.r