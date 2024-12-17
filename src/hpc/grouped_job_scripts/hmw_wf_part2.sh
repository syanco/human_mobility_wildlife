#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --job-name hmw_wf_part2
#SBATCH -c 4
##SBATCH --mem 500GB  
##SBATCH --mem-per-cpu 100G
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu

# set up paths
export wd=/scratch/julietcohen/covid_movement/human_mobility_wildlife
export src=$wd/src/workflow
cd $wd

# copy database to /tmp on worker node
cp $wd/processed_data/intermediate_db_copies/mosey_mod_clean-movement_complete.db /tmp/mosey_mod.db

# load R and spatial modules 
module load R/4.1.3 gdal/2.2.3 proj/5.2


# ------------ HPC step 8: Calculate Responses - Space Use - Fit dBBMMs ------------

# Make log file to track successful outputs
echo "species, ind_id, study_id, year, out_type, filename, produced, out_date" > $wd/out/dbbmm_log.csv
  
# Make big mem log file to track ind-year combos saved for the big mem parition
echo "species, ind_id, study_id, year, out_type, filename, produced, out_date" > $wd/out/dbbmm_bigmem_log.csv

Rscript $src/fit-dBBMMs.r /tmp/mosey_mod.db $wd/out 2

echo "fit-dBBMMs.r complete. Archiving database in /scratch."


# ------------ HPC step 9: Calculate Responses - Space Use - Calculate dBBMM areas and collate environment ------------

echo "species, ind_id, study_id, year, wk, area, sg, ghm, cbg_area, ndvi, tmax, n, a_bb, fixmed, m_error" > $wd/out/dbbmm_size.csv

Rscript $src/calc-space-use.r $wd/out /tmp/mosey_mod.db $wd/out/dbbmm_log.csv 4 -c

echo "calc-space-use.r complete."

# ------------ HPC step 10: Calculate Responses - Niche Breadth - Calculate MVNH Breadth ------------

# Create csv to store results with column names specified beforehand
# (this will overwrite the existing CSV if it already exists)
echo "total, tmax, ndvi, elev, cor, week, individual, scientificname, studyid, year, n" > $wd/out/niche_determinant_anthropause.csv
# Make log file to track successful outputs
echo "studyid, individual, scientificname, year, status, week" > $wd/out/niche_log.csv

# Execute calc size script
Rscript $src/calc-niche-breadth.r /tmp/mosey_mod.db ./out/niche_determinant_anthropause.csv 4

echo "calc-niche-breadth.r complete."

# ------------ HPC step 11: Inferential Models - Fit space use models part 1 ------------

# the skeleton for output dbbmm_size.csv was created prior to a previous step (Calculate dBBMM areas and collate environment)
Rscript $src/fit-space-use-dot-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area_dot 4 5 10000 5

echo "fit-space-use-dot-models.r complete."

# ------------ HPC step 12: Inferential Models - Fit space use models part 2 ------------

Rscript $src/fit-space-use-interactive-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area_interactive 4 5 10000 5

echo "fit-space-use-interactive-models.r complete."

# ------------ HPC step 13: Inferential Models - Fit space use models part 3 ------------

Rscript $src/fit-space-use-additive-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area_additive 4 5 10000 5

echo "fit-space-use-additive-models.r complete."

echo "SCRIPT COMPLETE: hmw_sf_part2"


