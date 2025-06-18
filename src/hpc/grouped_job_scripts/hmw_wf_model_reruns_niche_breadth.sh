#!/bin/bash

#SBATCH -t 03:00:00
#SBATCH --job-name reruns_niche_breadth
#SBATCH -c 3
#SBATCH --mem 200GB  
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu

# set up paths
export wd=/home/julietcohen/repositories/human_mobility_wildlife
export src=$wd/src/workflow
cd $wd

# copy database to /tmp on worker node
cp $wd/processed_data/intermediate_db_copies/mosey_mod_clean-movement_complete.db /tmp/mosey_mod.db

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

Rscript $src/refit-niche-breadth-additive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models_reruns/niche_additive 1 

echo "SCRIPT COMPLETE: refit-niche-breadth-additive-models.r"


# ------------ Model reruns for problematic MCMCs - Interactive Models ------------

echo "STARTING SCRIPT: refit-niche-breadth-interactive-models.r"

Rscript $src/refit-niche-breadth-interactive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models_reruns/niche_interactive 1 

echo "SCRIPT COMPLETE: refit-niche-breadth-interactive-models.r"

# ------------ Model summary for reruns ------------

echo "STARTING SCRIPT: niche_model_summaries_reruns.r"

Rscript $src/niche_model_summaries_reruns.r

echo "SCRIPT COMPLETE: niche_model_summaries_reruns.r"

echo "JOB COMPLETE: reruns_niche_breadth"