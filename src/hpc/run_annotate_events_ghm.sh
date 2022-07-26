#!/bin/bash

#SBATCH --job-name=ghm_annotation
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scott.yanco@yale.edu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100g -t 2:00:00
#SBATCH --partition=day

# module load R
module load miniconda
conda activate covid

# Declare WD
wd=~/project/covid-19_movement
src=$wd/analysis/src/workflow

# Move to WD
cd $wd

Rscript $src/annotate-events-ghm.r