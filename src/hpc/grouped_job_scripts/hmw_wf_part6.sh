#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --job-name hmw_wf_part6
#SBATCH -c 8
##SBATCH --mem 500GB  
##SBATCH --mem-per-cpu 100G
#SBATCH --mail-type ALL
#SBATCH --mail-user jscohen@ucsb.edu


# set up paths
export wd=/scratch/julietcohen/covid_movement/human_mobility_wildlife
export src=$wd/src/workflow
cd $wd

# load R and spatial modules 
module load R/4.1.3 gdal/2.2.3 proj/5.2

# ------------ HPC step 21: Intra-Individual Analysis - Area part 1 ------------

Rscript $src/fit_intra_ind_mod_additive_space.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 4 10000 5

echo "fit_intra_ind_mod_additive_space.r complete."


# ------------ HPC step 22: Intra-Individual Analysis - Area part 2 ------------

Rscript $src/fit_intra_ind_mod_interactive_space.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 4 10000 5

echo "fit_intra_ind_mod_interactive_space.r complete."


# ------------ HPC step 23: Intra-Individual Analysis - Niche part 1 ------------

Rscript $src/fit_intra_ind_mod_additive_niche.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 4 10000 5

echo "fit_intra_ind_mod_additive_niche.r complete."


# ------------ HPC step 24: Intra-Individual Analysis - Niche part 2 ------------

Rscript $src/fit_intra_ind_mod_interactive_niche.r $wd/out/dbbmm_size.csv $wd/out/intra_ind_models 4 10000 5

echo "fit_intra_ind_mod_interactive_niche.r complete."

# scp -r grace:$wdr/out/intra_ind_models $wd/out

echo "SCRIPT COMPLETE: hmw_sf_part6"