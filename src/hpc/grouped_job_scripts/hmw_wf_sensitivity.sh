#!/bin/bash

#SBATCH -t 48:00:00
#SBATCH --job-name hmw_wf_subsamples
##SBATCH -c 32
#SBATCH --mem 300G
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

# ------------ Niche breadth subsample of 50 ------------

echo "STARTING SCRIPT: calc-niche-breadth-subsample.r for sample size 50"

echo "total, tmax, ndvi, elev, cor, week, individual, scientificname, studyid, year, n" > $wd/out/niche_subsamples/niche_determinant_anthropause_subsample_50.csv

echo "studyid, individual, scientificname, year, status, week" > $wd/out/niche_subsamples/niche_subsample_log_50.csv

Rscript $src/calc-niche-breadth-subsample.r /tmp/mosey_mod.db ./out/niche_subsamples/niche_determinant_anthropause_subsample_50.csv 24 50

echo "SCRIPT COMPLETE: calc-niche-breadth-subsample.r for sample size 50"

# ------------ Niche breadth subsample of 40 ------------

echo "STARTING SCRIPT: calc-niche-breadth-subsample.r for sample size 40"

echo "total, tmax, ndvi, elev, cor, week, individual, scientificname, studyid, year, n" > $wd/out/niche_subsamples/niche_determinant_anthropause_subsample_40.csv

echo "studyid, individual, scientificname, year, status, week" > $wd/out/niche_subsamples/niche_subsample_log_40.csv

Rscript $src/calc-niche-breadth-subsample.r /tmp/mosey_mod.db ./out/niche_subsamples/niche_determinant_anthropause_subsample_40.csv 24 40

echo "SCRIPT COMPLETE: calc-niche-breadth-subsample.r for sample size 40"

# ------------ Niche breadth subsample of 30 ------------

echo "STARTING SCRIPT: calc-niche-breadth-subsample.r for sample size 30"

echo "total, tmax, ndvi, elev, cor, week, individual, scientificname, studyid, year, n" > $wd/out/niche_subsamples/niche_determinant_anthropause_subsample_30.csv

echo "studyid, individual, scientificname, year, status, week" > $wd/out/niche_subsamples/niche_subsample_log_30.csv

Rscript $src/calc-niche-breadth-subsample.r /tmp/mosey_mod.db ./out/niche_subsamples/niche_determinant_anthropause_subsample_30.csv 24 30

echo "SCRIPT COMPLETE: calc-niche-breadth-subsample.r for sample size 30"

# ------------ Niche breadth subsample of 20 ------------

echo "STARTING SCRIPT: calc-niche-breadth-subsample.r for sample size 20"

echo "total, tmax, ndvi, elev, cor, week, individual, scientificname, studyid, year, n" > $wd/out/niche_subsamples/niche_determinant_anthropause_subsample_20.csv

echo "studyid, individual, scientificname, year, status, week" > $wd/out/niche_subsamples/niche_subsample_log_20.csv

Rscript $src/calc-niche-breadth-subsample.r /tmp/mosey_mod.db ./out/niche_subsamples/niche_determinant_anthropause_subsample_20.csv 24 20

echo "SCRIPT COMPLETE: calc-niche-breadth-subsample.r for sample size 20"

# ------------ Niche breadth subsample of 10 ------------

echo "STARTING SCRIPT: calc-niche-breadth-subsample.r for sample size 10"

echo "total, tmax, ndvi, elev, cor, week, individual, scientificname, studyid, year, n" > $wd/out/niche_subsamples/niche_determinant_anthropause_subsample_10.csv

echo "studyid, individual, scientificname, year, status, week" > $wd/out/niche_subsamples/niche_subsample_log_10.csv

Rscript $src/calc-niche-breadth-subsample.r /tmp/mosey_mod.db ./out/niche_subsamples/niche_determinant_anthropause_subsample_10.csv 24 10

echo "SCRIPT COMPLETE: calc-niche-breadth-subsample.r for sample size 10"

# ------------ Plot niche breadth subsample ------------

# TODO: ADJUST FOLLOWING TO PLOT ALL 5 ITERATIONS OF SUBSAMPLES

#echo "STARTING SCRIPT: graph_niche_subsample.R"

#Rscript $src/graph_niche_subsample.R

#echo "SCRIPT COMPLETE: graph_niche_subsample.R"