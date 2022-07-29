#!/bin/bash

#SBATCH -t 12:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition day
#SBATCH -c 10
#SBATCH --mem-per-cpu 10G
#SBATCH -J calc_niches

# Load conda env
module load miniconda
conda activate covid

# Declare WD
wd=~/project/covid-19_movement
src=$wd/analysis/src/workflow

# Move to WD
cd $wd

#copy db to tmp
cp $wd/processed_data/mosey_mod.db /tmp/


# Execute calc size script/
Rscript $src/calc-niche-breadth.r /tmp/mosey_mod.db calc-niche-breadth 24

