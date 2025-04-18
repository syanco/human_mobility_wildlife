#!/bin/bash

#SBATCH -t 10:00:00
#SBATCH --job-name hmw_wf_part1
#SBATCH --mem 300G
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu


# set up paths
export wd=/home/julietcohen/repositories/human_mobility_wildlife
export src=$wd/src/workflow
cd $wd

# copy database to /tmp on worker node
cp $wd/processed_data/mosey_mod.db /tmp/mosey_mod.db

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


# ------------ HPC step 1: Filtering for Minimums ------------

echo "STARTING SCRIPT: filter_data_mins.R" 

Rscript $src/filter_data_mins.R --db /tmp/mosey_mod.db 30 3

echo "SCRIPT COMPLETE: filter_data_mins.R. Archiving database in /scratch."

# copy intermediate version of modified database to perisistent storage
cp /tmp/mosey_mod.db $wd/processed_data/intermediate_db_copies/mosey_mod_filter-data-mins_complete.db


# ------------ HPC step 2: Intersect events with Census Block Group geometries ------------

echo "STARTING SCRIPT: intersect-events-cbg.r"

Rscript $src/intersect-events-cbg.r

echo "SCRIPT COMPLETE: intersect-events-cbg.r."


# ------------ HPC step 3: Compute Census Block Group area ------------

echo "STARTING SCRIPT: compute-cbg-area.r"

Rscript $src/compute-cbg-area.r

echo "SCRIPT COMPLETE: compute-cbg-area.r."


# ------------ HPC step 4: Annotate Census Block Groups ------------

echo "STARTING SCRIPT: annotate-events-cbg.r"

Rscript $src/annotate-events-cbg.r

echo "SCRIPT COMPLETE: annotate-events-cbg.r."


# ------------ HPC step 5: Annotate events with SafeGraph ------------

echo "STARTING SCRIPT: annotate-events-safegraph.r"

Rscript $src/annotate-events-safegraph.r

echo "SCRIPT COMPLETE: annotate-events-safegraph.r."


# ------------ HPC step 6: Annotate Events with Human Modification ------------

echo "STARTING SCRIPT: annotate-events-ghm.r"

Rscript $src/annotate-events-ghm.r

echo "SCRIPT COMPLETE: annotate-events-ghm.r."


# ------------ HPC step 7: Clean Data ------------

echo "STARTING SCRIPT: clean_movement.r"

Rscript $src/clean_movement.r --db /tmp/mosey_mod.db

echo "SCRIPT COMPLETE: clean_movement.r. Archiving database in /scratch."

cp /tmp/mosey_mod.db $wd/processed_data/intermediate_db_copies/mosey_mod_clean-movement_complete.db

echo "JOB COMPLETE: hmw_wf_part1"
