#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --job-name hmw_wf_part5
#SBATCH -c 10
##SBATCH --mem 500GB  
##SBATCH --mem-per-cpu 100G
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu

# NOTE:
# These scripts only need to be run if the models from the prior workflow steps
# need to be rerun, according to the check of outputs from niche_model_summaries.r


# set up paths
export wd=/scratch/julietcohen/covid_movement/human_mobility_wildlife
export src=$wd/src/workflow
cd $wd

# load R and spatial modules 
module load R/4.1.3 gdal/2.2.3 proj/5.2

# ------------ HPC step 19: Model reruns for problematic MCMCs - Additive Models ------------

Rscript $src/refit-niche-breadth-additive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_additive 8 

echo "refit-niche-breadth-additive-models.r complete."


# ------------ HPC step 20: Model reruns for problematic MCMCs - Interactive Models ------------

Rscript $src/refit-niche-breadth-interactive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_interactive 8 

echo "refit-niche-breadth-interactive-models.r complete."

echo "SCRIPT COMPLETE: hmw_sf_part5"