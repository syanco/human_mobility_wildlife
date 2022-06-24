#!/bin/bash

#SBATCH --job-name=cbg_area
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ruth.oliver@yale.edu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100g -t 2-
#SBATCH --partition=general,pi_jetz
#SBATCH -C avx2

module load R/4.1.0-foss-2020b

# Declare WD
wd=~/project/covid-19_movement
src=$wd/analysis/src/workflow

# Move to WD
cd $wd

Rscript $src/compute-cbg-area.r