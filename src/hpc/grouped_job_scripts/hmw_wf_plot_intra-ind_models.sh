#!/bin/bash

#SBATCH -t 24:00:00
#SBATCH --job-name hmw_wf_plot_intra-ind_models
#SBATCH -c 10
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


# ------------ Plot intra-ind models space use ------------

Rscript $src/plot_intra-ind_model_area.r

echo "plot_intra-ind_model_area.r complete."


# ------------ Plot intra-ind models niche breadth ------------

Rscript $src/plot_intra-ind_model_niche.r

echo "plot_intra-ind_model_niche.r complete."

echo "JOB COMPLETE: hmw_wf_plot_intra-ind_models"
