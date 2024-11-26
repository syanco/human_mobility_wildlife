#!/bin/bash

#SBATCH -t 2:00:00                    # max time for job
#SBATCH --mail-type ALL               # email all job events
#SBATCH --mail-user jscohen@ucsb.edu
##SBATCH --partition day              # partition "day" may only exist on Yale HPC
#SBATCH -c 24                         # number of CPU's for job
##SBATCH --mem-per-cpu 10G            # can add this back in if need to specify memory
#SBATCH -J calc_niches

# Declare WD
wd=/scratch/julietcohen/covid_movement/human_mobility_wildlife
src=$wd/src/workflow

cd $wd

#copy db to tmp
cp $wd/processed_data/mosey_mod.db /tmp/


# Execute calc size script
Rscript $src/calc-niche-breadth.r /tmp/mosey_mod.db ./out/niche_determinant_anthropause.csv 24

