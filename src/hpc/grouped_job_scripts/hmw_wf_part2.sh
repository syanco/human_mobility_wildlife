#!/bin/bash

#SBATCH -t 96:00:00
#SBATCH --job-name hmw_wf_part2
#SBATCH -c 32
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu
#SBATCH --nodelist=node53
#SBATCH --partition=largemem
#SBATCH --mem=900G 

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


# ------------ HPC step 8: Calculate Responses - Space Use - Fit dBBMMs ------------

# Make log file to track successful outputs
echo "species, ind_id, study_id, year, out_type, filename, produced, out_date" > $wd/out/dbbmm_log.csv
  
# Make big mem log file to track ind-year combos saved for the big mem parition
echo "species, ind_id, study_id, year, out_type, filename, produced, out_date" > $wd/out/dbbmm_bigmem_log.csv

# Make log file to track how many individual and year pairs do not exist in the event data
echo "study_id, individual_id, year" > $wd/out/no_ind_yr_pairs.csv

echo "STARTING SCRIPT: fit-dBBMMs.r"

Rscript $src/fit-dBBMMs.r /tmp/mosey_mod.db $wd/out 7

echo "SCRIPT COMPLETE: fit-dBBMMs.r"


# ------------ HPC step 9: Calculate Responses - Space Use - Calculate dBBMM areas and collate environment ------------

echo "species, ind_id, study_id, year, wk, area, sg, ghm, cbg_area, ndvi, tmax, n, a_bb, fixmed, m_error" > $wd/out/dbbmm_size.csv

echo "STARTING SCRIPT: calc-space-use.r"

Rscript $src/calc-space-use.r $wd/out /tmp/mosey_mod.db $wd/out/dbbmm_log.csv 7 -c

echo "SCRIPT COMPLETE: calc-space-use.r"


# ------------ HPC step 10: Calculate Responses - Niche Breadth - Calculate MVNH Breadth ------------

# Create csv to store results with column names specified beforehand
# (this will overwrite the existing CSV if it already exists)
echo "total, tmax, ndvi, elev, cor, week, individual, scientificname, studyid, year, n" > $wd/out/niche_determinant_anthropause.csv
# Make log file to track successful outputs
echo "studyid, individual, scientificname, year, status, week" > $wd/out/niche_log.csv

echo "STARTING SCRIPT: calc-niche-breadth.r"

# Execute calc size script
Rscript $src/calc-niche-breadth.r /tmp/mosey_mod.db ./out/niche_determinant_anthropause.csv 7

echo "SCRIPT COMPLETE: calc-niche-breadth.r"


# ------------ HPC step 11: Inferential Models - Fit space use models part 1 ------------

echo "STARTING SCRIPT: fit-space-use-dot-models.r"

Rscript $src/fit-space-use-dot-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area_dot 7 3 10000 5

echo "SCRIPT COMPLETE: fit-space-use-dot-models.r"


# ------------ HPC step 12: Inferential Models - Fit space use models part 2 ------------

echo "STARTING SCRIPT: fit-space-use-interactive-models.r"

Rscript $src/fit-space-use-interactive-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area_interactive 7 3 10000 5

echo "SCRIPT COMPLETE: fit-space-use-interactive-models.r"


# ------------ HPC step 13: Inferential Models - Fit space use models part 3 ------------

echo "STARTING SCRIPT: fit-space-use-additive-models.r"

Rscript $src/fit-space-use-additive-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area_additive 7 3 10000 5

echo "SCRIPT COMPLETE: fit-space-use-additive-models.r"


# ------------ HPC step 14: Visualize model outputs from space use models ------------

echo "STARTING SCRIPT: area_model_summaries.r"

Rscript $src/area_model_summaries.r

echo "SCRIPT COMPLETE: area_model_summaries.r."

echo "JOB COMPLETE: hmw_sf_part2"