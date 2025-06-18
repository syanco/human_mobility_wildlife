#!/bin/bash

#SBATCH -t 10:00:00
#SBATCH --job-name fix_rate
#SBATCH -c 10
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu
#SBATCH --mem=250G 

# set up paths
export wd=/home/julietcohen/repositories/human_mobility_wildlife
export src=$wd/src/workflow
cd $wd

# copy database to /tmp on worker node
cp $wd/processed_data/intermediate_db_copies/mosey_mod_clean-movement_complete.db /tmp/mosey_mod.db

module load R/4.3.1

# set paths for PROJ and GDAL installations in my custom conda env
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

echo "study_id, species, individual_id, year, fix_rate_med, fix_rate_min, fix_rate_max" > $wd/out/fixrate_med_min_max.csv

Rscript $src/fix_rate_fig.r /tmp/mosey_mod.db $wd/out 7

echo "JOB COMPLETE"