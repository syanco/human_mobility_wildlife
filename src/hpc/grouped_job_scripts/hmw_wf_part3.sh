#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --job-name hmw_wf_part3
#SBATCH -c 5
##SBATCH --mem 500GB  
##SBATCH --mem-per-cpu 100G
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu

# NOTE:
# These scripts only need to be run if the models from the prior workflow steps
# need to be rerun, according to the check of outputs from area_model_summaries.r

# set up paths
export wd=/scratch/julietcohen/covid_movement/human_mobility_wildlife
export src=$wd/src/workflow
cd $wd

# load R and spatial modules 
module load R/4.1.3 gdal/2.2.3 proj/5.2

# ------------ HPC step 14: Model reruns for problematic MCMCs - Additive Models ------------

Rscript $src/refit-space-use-additive-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area_additive 4

echo "refit-space-use-additive-models.r complete."

# ------------ HPC step 15: Model reruns for problematic MCMCs - Interactive Models ------------

Rscript $src/refit-space-use-interactive-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area_additive 4

echo "refit-space-use-interactive-models.r complete."

echo "SCRIPT COMPLETE: hmw_sf_part3"


