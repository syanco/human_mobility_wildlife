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


# ------------ Plot intra-ind models space use ------------

Rscript $src/plot_intra-ind_model_area.r

echo "plot_intra-ind_model_area.r complete."


# ------------ Plot intra-ind models niche breadth ------------

Rscript $src/plot_intra-ind_model_niche.r

echo "plot_intra-ind_model_niche.r complete."

echo "JOB COMPLETE: hmw_wf_plot_intra-ind_models"
