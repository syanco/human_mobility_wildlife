#!/bin/bash

#SBATCH -t 04:00:00
#SBATCH --job-name hmw_wf_pre_anno
#SBATCH --mem 300GB  
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu


# set up paths
export wd=/home/julietcohen/repositories/human_mobility_wildlife
export src=$wd/src/workflow
cd $wd

# copy database from mosey output to human_mobility_wildlife/raw_data
cp /home/julietcohen/repositories/mosey_db_output/data/mosey.db $wd/raw_data/mosey.db
# copy database to /tmp on worker node
cp $wd/raw_data/mosey.db /tmp/mosey_mod.db

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


# ------------ Remove some skunk data ------------

echo "STARTING SCRIPT: adjust_skunk_data.R"

Rscript $src/adjust_skunk_data.R

echo "SCRIPT COMPLETE: adjust_skunk_data.R"

# ------------ Trimming ------------

echo "STARTING SCRIPT: trim_data.r" 

Rscript $src/trim_data.r

echo "SCRIPT COMPLETE: trim_data.r"

cp /tmp/mosey_mod.db $wd/processed_data/mosey_mod.db

echo "JOB COMPLETE: hmw_wf_pre_anno"
