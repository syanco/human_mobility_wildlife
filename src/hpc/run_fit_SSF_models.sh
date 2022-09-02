#!/bin/bash

#SBATCH -t 1-
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition day
#SBATCH -c 24
#SBATCH --mem-per-cpu 10G
#SBATCH -J fit_SSFs

# Load conda env
module load miniconda
conda activate amt_db

# Declare WD
wd=/gpfs/loomis/pi/jetz/sy522/covid-19_movement

# Move to WD
cd $wd

#copy db to tmp
cp $wd/processed_data/mosey_mod.db /tmp/

# Execute cleaning script (non parallel)
Rscript $wd/analysis/src/workflow/fit-SSF-models.r /tmp/mosey_mod.db $wd/out/ssf-background-pts/annotated $wd/out/ssf-mods 24 10 10000 4
