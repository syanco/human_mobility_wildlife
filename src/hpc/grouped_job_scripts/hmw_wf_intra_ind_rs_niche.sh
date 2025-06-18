#!/bin/bash

#SBATCH -t 10:00:00
#SBATCH --job-name rs_n
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

# ------------ Niche Intra-Individual Interactive Analysis Random Slopes ------------

echo "STARTING SCRIPT: fit_intra_ind_int_mod_random_slopes_niche.r"

Rscript $src/fit_intra_ind_int_mod_random_slopes_niche.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 12 10000 5

echo "SCRIPT COMPLETE: fit_intra_ind_int_mod_random_slopes_niche.r"


# ------------ Niche Intra-Individual Additive Analysis Random Slopes ------------

echo "STARTING SCRIPT: fit_intra_ind_add_mod_random_slopes_niche.r"

Rscript $src/fit_intra_ind_add_mod_random_slopes_niche.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 12 10000 5

echo "SCRIPT COMPLETE: fit_intra_ind_add_mod_random_slopes_niche.r"

echo "JOB COMPLETE: rs_n"