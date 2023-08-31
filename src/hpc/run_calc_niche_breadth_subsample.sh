#!/bin/bash

#SBATCH -t 6:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user scott.yanco@yale.edu
#SBATCH --partition day
#SBATCH -c 24
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


# Vector of subsample sizes
number_vector=(40 30 20 10)  # Add your desired numbers here

# Loop over the vector and execute calc size script
for num in "${number_vector[@]}"; do
  # Execute calc size script/
  Rscript $src/calc-niche-breadth.r /tmp/mosey_mod.db ./out/niche_determinant_anthropause.csv 24 $v
done
