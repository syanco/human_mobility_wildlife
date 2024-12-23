#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --job-name hmw_wf_part3
#SBATCH -c 27
##SBATCH --mem 500GB  
##SBATCH --mem-per-cpu 100G
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu


# set up paths
export wd=/scratch/julietcohen/covid_movement/human_mobility_wildlife
export src=$wd/src/workflow
cd $wd

# load R and spatial modules 
module load R/4.1.3 gdal/2.2.3 proj/5.2


# ------------ HPC step 15: Fit niche breadth models part 1 ------------

echo "STARTING SCRIPT: fit-niche-breadth-dot-models.r"

Rscript $src/fit-niche-breadth-dot-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_dot 24 5 10000 5

echo "SCRIPT COMPLETE: fit-niche-breadth-dot-models.r"


# ------------ HPC step 16: Fit niche breadth models part 2 ------------

echo "STARTING SCRIPT: fit-niche-breadth-additive-models.r"

Rscript $src/fit-niche-breadth-additive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_additive 24 5 10000 5

echo "SCRIPT COMPLETE: fit-niche-breadth-additive-models.r"


# ------------ HPC step 17: Fit niche breadth models part 3 ------------

echo "STARTING SCRIPT: fit-niche-breadth-interactive-models.r"

Rscript $src/fit-niche-breadth-interactive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_interactive 24 5 10000 5

echo "SCRIPT COMPLETE: fit-niche-breadth-interactive-models.r"

# ------------ HPC step 18: Visualize model outputs from niche breadth models ------------

echo "STARTING SCRIPT: niche_model_summaries.r"

Rscript $src/niche_model_summaries.r

echo "SCRIPT COMPLETE: niche_model_summaries.r"

echo "JOB COMPLETE: hmw_sf_part3"