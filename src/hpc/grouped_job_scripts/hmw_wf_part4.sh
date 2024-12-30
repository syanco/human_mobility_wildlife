#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --job-name hmw_wf_part4
#SBATCH -c 8
##SBATCH --mem 500GB  
##SBATCH --mem-per-cpu 100G
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu


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


# ------------ HPC step 19: Intra-Individual Analysis - Area part 1 ------------

echo "STARTING SCRIPT: fit_intra_ind_mod_additive_space.r"

Rscript $src/fit_intra_ind_mod_additive_space.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 4 10000 5

echo "SCRIPT COMPLETE: fit_intra_ind_mod_additive_space.r"


# ------------ HPC step 20: Intra-Individual Analysis - Area part 2 ------------

echo "STARTING SCRIPT: fit_intra_ind_mod_interactive_space.r"

Rscript $src/fit_intra_ind_mod_interactive_space.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 4 10000 5

echo "SCRIPT COMPLETE: fit_intra_ind_mod_interactive_space.r"


# ------------ HPC step 21: Intra-Individual Analysis - Niche part 1 ------------

echo "STARTING SCRIPT: fit_intra_ind_mod_additive_niche.r"

Rscript $src/fit_intra_ind_mod_additive_niche.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 4 10000 5

echo "SCRIPT COMPLETE: fit_intra_ind_mod_additive_niche.r"


# ------------ HPC step 22: Intra-Individual Analysis - Niche part 2 ------------

echo "STARTING SCRIPT: fit_intra_ind_mod_interactive_niche.r"

Rscript $src/fit_intra_ind_mod_interactive_niche.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 4 10000 5

echo "SCRIPT COMPLETE: fit_intra_ind_mod_interactive_niche.r"

echo "JOB COMPLETE: hmw_sf_part4"