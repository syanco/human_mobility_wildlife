#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --job-name hmw_wf_part1
#SBATCH -c 27
#SBATCH --mem 500GB  
##SBATCH --mem-per-cpu 100G
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu

# set up paths
export wd=/scratch/julietcohen/covid_movement/human_mobility_wildlife
export src=$wd/src/workflow
cd $wd

# copy database to /tmp on worker node
cp $wd/processed_data/mosey_mod.db /tmp/

# load R and spatial modules 
module load R/4.1.3 gdal/2.2.3 proj/5.2


# ------------ HPC step 1: Filtering for Minimums ------------

Rscript $src/filter_data_mins.R --db /tmp/mosey_mod.db 30 5

echo "filer_data_mins.R complete. Archiving database in /scratch."

# copy intermediate version of modified database to perisistent storage
cp /tmp/mosey_mod.db $wd/processed_data/intermediate_db_copies/mosey_mod_filter-data-mins_complete.db


# ------------ HPC step 2: Intersect events with Census Block Group geometries ------------

Rscript $src/intersect-events-cbg.r

echo "intersect-events-cbg.r complete."


# ------------ HPC step 3: Compute Census Block Group area ------------

Rscript $src/compute-cbg-area.r

echo "compute-cbg-area.r complete."


# ------------ HPC step 4: Annotate Census Block Groups ------------

Rscript $src/annotate-events-cbg.r

echo "annotate-events-cbg.r complete. Archiving database in /scratch."

cp /tmp/mosey_mod.db $wd/processed_data/intermediate_db_copies/mosey_mod_annotate-events-cbg_complete.db


# ------------ HPC step 5: Annotate events with SafeGraph ------------

Rscript $src/annotate-events-safegraph.r

echo "annotate-events-safegraph.r complete. Archiving database in /scratch."

cp /tmp/mosey_mod.db $wd/processed_data/intermediate_db_copies/mosey_mod_annotate-events-safegraph_complete.db


# ------------ HPC step 6: Annotate Events with Human Modification ------------

Rscript $src/annotate-events-ghm.r

echo "annotate-events-ghm.r complete. Archiving database in /scratch."

cp /tmp/mosey_mod.db $wd/processed_data/intermediate_db_copies/mosey_mod_annotate-events-ghm_complete.db


# ------------ HPC step 7: Clean Data ------------

Rscript $src/clean_movement.r --db /tmp/mosey_mod.db

echo "clean_movement.r complete. Archiving database in /scratch."

cp /tmp/mosey_mod.db $wd/processed_data/intermediate_db_copies/mosey_mod_clean-movement_complete.db

echo "SCRIPT COMPLETE: hmw_sf_part1"
