#!/bin/bash

#SBATCH --job-name=safegraph_annotation
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ruth.oliver@yale.edu
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000g -t 2-
#SBATCH --partition=general,bigmem,pi_jetz
#SBATCH -C avx2

module load R

# Declare WD
wd=~/project/covid-19_movement
src=$wd/analysis/src/workflow

# Move to WD
cd $wd

Rscript $src/annotate-events-safegraph.r