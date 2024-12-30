#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --job-name reruns_niche_breadth
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

module load R/4.3.1

# set paths for PROJ, GDAL, and compiler installations in my conda env
source /home/julietcohen/miniconda3/etc/profile.d/conda.sh
conda activate r_spatial2
export LD_LIBRARY_PATH=/home/julietcohen/miniconda3/envs/r_spatial2/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/home/julietcohen/miniconda3/envs/r_spatial2/lib/pkgconfig:$PKG_CONFIG_PATH
export CXX=/home/julietcohen/miniconda3/envs/r_spatial2/bin/g++
export CC=/home/julietcohen/miniconda3/envs/r_spatial2/bin/gcc

# add logging to ensure the new shell is using the intended software
echo "GDAL version:"
gdal-config --version
echo "PROJ version:"
proj
echo "C++ compiler version:"
$CXX --version
echo "Make version:"
make --version


# ------------ Model reruns for problematic MCMCs - Additive Models ------------

echo "STARTING SCRIPT: refit-niche-breadth-additive-models.r"

Rscript $src/refit-niche-breadth-additive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_additive 8 

echo "SCRIPT COMPLETE: refit-niche-breadth-additive-models.r"


# ------------ Model reruns for problematic MCMCs - Interactive Models ------------

echo "STARTING SCRIPT: refit-niche-breadth-interactive-models.r"

Rscript $src/refit-niche-breadth-interactive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_interactive 8 

echo "SCRIPT COMPLETE: refit-niche-breadth-interactive-models.r"

echo "JOB COMPLETE: reruns_niche_breadth"