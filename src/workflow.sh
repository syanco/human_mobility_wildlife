                  ##########################################                    
                  ####                                  ####
                  ####    COVID-19 & Animal Movement    ####
                  ####          Workflow Script         ####
                  ####                                  ####
                  ####              contact:            ####
                  ####       scott.yanco@yale.edu       ####
                  ####                                  ####
                  ##########################################
                  
#------------------------------------------------------------------------------#

# Top-level workfow script to execute/describe steps to run analyses associated 
# with the North American COVID-19 & Animal Movemnt project.  Some steps are 
# only described (i.e., do not contain code) due to restrictions on data, access
# credentials, non-coded steps, etc.

# CONDA ENVIRONMENTS
# 
# Much of this code is meant to be executed within coda environments. The 
# directory /analysis/src/conda_envs contains .yml files from which the 
# requisite conda environments can be built.

#------------------------------------------------------------------------------#


####----  Initialization  ----####

# Define working directory - should point to project root folder
wd=~/projects/covid-19_movement


####----  Build Initial Database  ----####

# Go to working directory
cd $wd

# Build a mosey_db-style SQLite database from movebank to store raw data
# See https://github.com/benscarlson/mosey_db for more detail.

# TODO:  @diego & @ben, add 'final' code to build db here


####----  Prep and Clean Data ----####

# Make a copy of db to be modified
cp $wd/raw_data/mosey.db $wd/processed_data/mosey_mod.db

#activate spatial env
conda activate spatial

# Process and clean data:
#   * add study period annotation
#   * remove events outside study area
#   * basic data cleaning
Rscript $wd/analysis/src/01-prep_and_clean.r --db $wd/processed_data/mosey_mod.db
