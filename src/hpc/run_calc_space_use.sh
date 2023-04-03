#!/bin/bash

#SBATCH -t 22:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition day
#SBATCH -c 10
#SBATCH --mem-per-cpu 20G
#SBATCH -J calc_size_with_TMAX-TMIN

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
Rscript $src/calc-space-use.r ./out /tmp/mosey_mod.db ./out/dbbmm_log.csv 10 -c

