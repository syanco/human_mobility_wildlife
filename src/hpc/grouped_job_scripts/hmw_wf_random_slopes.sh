#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --job-name random_slopes
#SBATCH -c 15
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu
#SBATCH --mem=300G


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


# ------------ Intra-Individual Analysis Random Slopes Niche ------------

echo "STARTING SCRIPT: fit_intra_ind_mod_random_slopes_niche.r"

Rscript $src/fit_intra_ind_mod_random_slopes_niche.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 10 10000 5

echo "SCRIPT COMPLETE: fit_intra_ind_mod_random_slopes_niche.r"

# ------------ Intra-Individual Analysis Raondom Slopes Area ------------

echo "STARTING SCRIPT: fit_intra_ind_mod_random_slopes_space.r"

Rscript $src/fit_intra_ind_mod_random_slopes_space.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 10 20000 2

echo "SCRIPT COMPLETE: fit_intra_ind_mod_random_slopes_space.r"

echo "JOB COMPLETE: random_slopes"

# ------------ Intra-Individual Analysis - Niche rerun ------------

echo "STARTING SCRIPT: fit_intra_ind_mod_additive_niche.r"

Rscript $src/fit_intra_ind_mod_additive_niche.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 10 10000 5

echo "SCRIPT COMPLETE: fit_intra_ind_mod_additive_niche.r"

