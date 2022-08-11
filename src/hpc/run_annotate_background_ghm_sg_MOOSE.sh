#!/bin/bash

#SBATCH --job-name=ghm_annotation-moose-BG
#SBATCH --cpus-per-task=15
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scott.yanco@yale.edu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100g -t 3:00:00
#SBATCH --partition=bigmem

# module load R
module load miniconda
conda activate spatial

# Declare WD
wd=/gpfs/loomis/pi/jetz/sy522/covid-19_movement
src=$wd/analysis/src/workflow

# Move to WD
cd $wd

Rscript $src/annotate-background-sg-ghm-MOOSE_scratch.r 15