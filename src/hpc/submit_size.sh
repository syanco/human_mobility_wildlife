#!/bin/bash

#SBATCH -t 3:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition pi_jetz
#SBATCH -c 24
#SBATCH --mem-per-cpu 10G
#SBATCH -J calc_size_2022-02-11

# Load conda env
module load miniconda
conda activate covid

# Declare WD
wd=~/project/covid-19_movement

#copy db to tmp
cp $wd/processed_data/mosey_mod.db /tmp/


# Execute calc size script/
Rscript $wd/analysis/src/0X-calc_space_use.r /tmp/mosey_mod.db ./out/dbbmm_log.csv 24

