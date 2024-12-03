#!/bin/bash

#SBATCH -t 04:00:00                    # max time for job
#SBATCH --mail-type ALL               # email all job events
#SBATCH --mail-user jscohen@ucsb.edu
##SBATCH --partition day              # partition "day" may only exist on Yale HPC
#SBATCH -c 10                         # number of CPU's, increase later if needed, <20 for high mem nodes on Pod 
##SBATCH --mem-per-cpu 10G            # can add this back in if need to specify memory
#SBATCH -J calc_niches
#SBATCH --nodes 1

# Declare WD
wd=/scratch/julietcohen/covid_movement/human_mobility_wildlife
src=$wd/src/workflow

# Create csv to store results with column names specified beforehand
# (this will overwrite the existing CSV if it already exists)
echo "total, tmax, ndvi, elev, cor, week, individual, scientificname, studyid, year, n" > $wd/out/niche_determinant_anthropause.csv
# Make log file to track successful outputs
echo "studyid, individual, scientificname, year, status, week" > $wd/out/niche_log.csv

#copy db to tmp
cp $wd/processed_data/mosey_mod.db /tmp/

# Execute calc size script
Rscript $src/calc-niche-breadth.r /tmp/mosey_mod.db ./out/niche_determinant_anthropause.csv 10

