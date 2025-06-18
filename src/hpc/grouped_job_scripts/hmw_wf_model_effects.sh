#!/bin/bash

#SBATCH -t 00:30:00
#SBATCH --job-name hmw_wf_model_effects
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


# ------------ Space Use Model Effects ------------

echo "STARTING SCRIPT: select_space_use_model-effects.r" 

Rscript $src/select_space_use_model-effects.r $wd/out/single_species_models $wd/out/figs /home/julietcohen/covid_movement_full_repo/raw_data/anthropause_data_sheet.csv

echo "SCRIPT COMPLETE: select_space_use_model-effects.r"


# ------------ Area Effects ------------

echo "STARTING SCRIPT: estimate_area_effects.r" 

Rscript $src/estimate_area_effects.r $wd/out/single_species_models $wd/out/figs

echo "SCRIPT COMPLETE: estimate_area_effects.r"


# ------------ Niche Model Effects ------------

echo "STARTING SCRIPT: select_niche_model_effects.r" 

Rscript $src/select_niche_model-effects.r $wd/out/single_species_models $wd/out/figs /home/julietcohen/covid_movement_full_repo/raw_data/anthropause_data_sheet.csv

echo "SCRIPT COMPLETE: select_niche_model_effects.r"


# ------------ Niche Effects ------------

echo "STARTING SCRIPT: estimate_niche_effects.r" 

Rscript $src/estimate_niche_effects.r $wd/out/single_species_models $wd/out/figs

echo "SCRIPT COMPLETE: estimate_niche_effects.r"
