#!/bin/bash

#SBATCH -t 48:00:00
#SBATCH --job-name hmw_wf_part3
#SBATCH -c 1
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu
#SBATCH --mem=500G 
##SBATCH --nodelist=node53
##SBATCH --partition=largemem

# set up paths
export wd=/home/julietcohen/repositories/human_mobility_wildlife
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


# ------------ HPC step 15: Fit niche breadth models part 1 ------------

echo "STARTING SCRIPT: fit-niche-breadth-dot-models.r"

Rscript $src/fit-niche-breadth-dot-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_dot 7 3 10000 5

echo "SCRIPT COMPLETE: fit-niche-breadth-dot-models.r"


# ------------ HPC step 16: Fit niche breadth models part 2 ------------

echo "STARTING SCRIPT: fit-niche-breadth-additive-models.r"

Rscript $src/fit-niche-breadth-additive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_additive 7 3 10000 5

echo "SCRIPT COMPLETE: fit-niche-breadth-additive-models.r"


# ------------ HPC step 17: Fit niche breadth models part 3 ------------

echo "STARTING SCRIPT: fit-niche-breadth-interactive-models.r"

Rscript $src/fit-niche-breadth-interactive-models.r $wd/out/niche_determinant_anthropause.csv $wd/out/dbbmm_size.csv $wd/out/single_species_models/niche_interactive 7 3 10000 5

echo "SCRIPT COMPLETE: fit-niche-breadth-interactive-models.r"

# ------------ HPC step 18: Visualize model outputs from niche breadth models ------------

echo "STARTING SCRIPT: niche_model_summaries.r"

Rscript $src/niche_model_summaries.r

echo "SCRIPT COMPLETE: niche_model_summaries.r"

echo "JOB COMPLETE: hmw_sf_part3"