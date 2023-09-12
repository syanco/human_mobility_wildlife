                  ##########################################                    
                  ####                                  ####
                  ####    COVID-19 & Animal Movement    ####
                  ####          Workflow Script         ####
                  ####                                  ####
                  ####             contact:             ####
                  ####       scott.yanco@yale.edu       ####
                  ####                                  ####
                  ##########################################
                  
#------------------------------------------------------------------------------#

# Top-level workfow script to executeanalyses associated 
# with the North American COVID-19 & Animal Movement project. 

# This project is led by: Ruth Oliver & Scott Yanco

# This workflow is generally built to be run in a high performance computing 
# (HPC) environment. Thus, most steps in the workflow simply use `sbatch` to 
# submit a job to the slurm manager.  However, we also include (commented) code 
# to run each step without submission to a slurm manager if desired 
# (not recommended).

# This workflow also assumes a specific directory strcuture.  See README for 
# details on directory strcuture.

# CONDA ENVIRONMENTS
# 
# Much of this code is meant to be executed within coda environments. The 
# directory /analysis/src/conda_envs contains .yml files from which the 
# requisite conda environments can be built. See readme for details on building 
# conda environments

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#

  # TODO:

#------------------------------------------------------------------------------#

####----  Initialization  ----####
  
  # define local working directory
  wd=/Users/scottyanco/Documents/covid-19_movement
  src=$wd/analysis/src
  cd $wd
    
  # defiune remot (HPC working directory)
  wdr=/home/sy522/project/covid-19_movement
  srcr=$wdr/analysis/src
  
  # Go to working directory
  cd $wdr
  
####



####----  Prep and Clean Data ----####

  ##-- Working DB --##

    # Make a copy of db to be modified
    cp $wd/raw_data/mosey.db $wd/processed_data/mosey_mod_2023.db
    # cp $wd/raw_data/mosey_swap.db $wd/processed_data/mosey_swap_mod.db
  
  ##
  
  ##-- Trim Data to Time/Area of Interest --##
  
  Rscript $src/workflow/trim_data.r --db $wd/processed_data/mosey_mod_2023.db


####----  Annotate Data ----####

  ##-- Environmental Annotations --#
  
    # Inputs: db:event_clean + analysis/ctfs/env.csv
    # Outputs: db:event_clean (annotated)
    
    # INTERACTIVE - RUN LOCAL (script must be run interactively; pulls db from hpc then puts it back)
    $src/workflow/wf-mosey_env.sh
  
  ##
    
  #**** Move to HPC here. ****#
  
  ##-- Census Annotations --##
  
    #- Intersect events with cbg geometries -#
    
      # Inputs:   db:event_clean + cbg shp file
      # Outputs:  csv per job (event_id + cbg info)

      # DSQ: (Yale-specific)
      # module load R/4.1.0-foss-2020b
      conda activate covid
      module load dSQ
      
      Rscript $srcr/workflow/create_intersection_joblist.r

      dsq --job-file $srcr/workflow/joblist.txt --mem-per-cpu 40g -t 02:00:00

      # The step above generates a .sh file to submit the job to the Slurm manager
      # Thus, after running the previous line, thje file referenced below will be 
      # created (and update the date to match the day it was generated).
      # sbatch dsq-joblist-2022-07-22.sh
      sbatch dsq-joblist-2023-09-08.sh
    
    #


    #- Compute cbg area -#

      # Inputs: cbg shp files
      # Outputs: csv (cbg info + area)

      # SLURM:
      sbatch $srcr/hpc/run_compute_cbg_area.sh
      
      # ON DEMAND:
      # Rscript $src/workflow/compute-cbg-area.r
      
    #
    
    #- Annotate events with cbg info -#
    
      # Inputs: db:event_clean + cbg intersection csv + cbg area csv
      # Outputs: csv (event_id + cbg info + cbg area)

      # SLURM:
      sbatch $srcr/hpc/run_annotate_events_cbg.sh
      
      # ON DEMAND:
      # Rscript $src/workflow/annotate-events-cbg.r
    
    #
  
  ##


  ##-- SafeGraph Annotations --##
  
    #- Process SafeGraph data -#
      
      # Inputs: safegraph txt files (one file per cbg/week)
      # Outputs: sg data csv (one file per cbg/week)

      # SLURM:
      sbatch $srcr/hpc/run_process_safegraph_data.sh
      
      # ON DEMAND:
      # Rscript $src/workflow/process-safegraph-data.r
    
    #

    #- Annotate events with SafeGraph -#
    
      # Inputs: db:event_clean + cbg info csv + sg data csv
      # Outputs: csv (event_id + timestamp + cbg info = cbg area + sg count)
  
      # SLURM:
      sbatch $srcr/hpc/run_annotate_events_safegraph.sh
      
      # ON DEMAND:
      # Rscript $src/workflow/annotate-events-safegraph.r
  
    #

  ##


  ##-- Human Modification Annotations --#
  
    #- Annotate events with TNC GHM -#
    
      # Inputs: db:event_clean + ghm raster
      # Ouputs: csv (event_id + ghm)

      # SLURM:
      sbatch $srcr/hpc/run_annotate_events_ghm.sh
    
      # ON DEMAND:
      # Rscript $src/workflow/annotate-events-ghm.r
    
    #
  
  ##

  ##-- Clean Data --##
    
      # Inputs:   db:event_trim 
      # Outputs:  db:event_mod + db:event_clean
  
    
      # SLURM:
      sbatch $srcr/hpc/run_clean_movement.sh
  
      # ON DEMAND:
      # conda activate spatial
      # Rscript $src/workflow/clean_movement.r --db $wd/processed_data/mosey_mod.db
    
  ##- MVNH Breadth Subsampling -#
      
      # Create csv to store results
      echo "total, tmax, ndvi, elev, cor, week, individual, scientificname, studyid, year, n" > ./out/niche_determinant_anthropause_subsample.csv

      # Make log file to track successful outputs
      echo "studyid, individual, scientificname, year, status, week" > out/niche_log_subsample.csv
  
      # Inputs: db:event_clean  + out filepath + no. cores
      # Outputs: csv 
      
      sbatch $srcr/hpc/run_calc_niche_breadth_subsample.sh

  ##

 ##- Filter to Min Sample - Report Sample Size  -#
      

      sbatch $srcr/hpc/run_filter_data_mins.sh

  ##



####

####----  Calculate Repsonses ----####

  ##---  Space Use  ---##

    #- Fit dBBMMs -#
      
      # Inputs: db:event_clean + out/dbbmm_log.csv (blank) + 
      #           out/dbbmm_bigmem_log.csv (blank) + no. cores
      # Outputs: rdata files in out/ (one per individual) 

      # Make log file to track successful outputs
      echo "species, ind_id, study_id, year, out_type, filename, produced, out_date" > out/dbbmm_log.csv
  
      # Make big mem log file to track ind-year combos saved for the big mem parition
      echo "species, ind_id, study_id, year, out_type, filename, produced, out_date" > out/dbbmm_bigmem_log.csv
  
      # SLURM:
      sbatch ~/project/covid-19_movement/analysis/src/hpc/run_fit_dBBMMs.sh

      # ON DEMAND:
      # conda activate covid
      # Rscript $src/workflow/fit-dBBMMs.r $wd/processed_data/mosey_mod.db $wd/out 1
      
      # TODO:
      #   * big mem log likely unneccesary, deprecate/delete
  
    #


    #- Calculate dBBMM areas and collate environment -#

      # Inputs: db:event_clean + db:event_census + out filepath + 
      #           out/dbbmm_log.csv (filled) + no. cores
      # Outputs: csv 
      
      # Create csv to store results
      echo "species, ind_id, study_id, year, wk, area, sg, ghm, cbg_area, ndvi, tmax, n, a_bb, fixmed, m_error" > ./out/dbbmm_size.csv

      # SLURM:
      sbatch $src/hpc/run_calc_space_use.sh
      
      # ON DEMAND:
      # conda activate covid
      # Rscript $src/workflow/calc-space-use.r
    
    #

    #- Check area size ~ sample size -#
    
    Rscript $src/workflow/check_area_size_sample_balance.R

    #

  ##
  
  
  ##-- Niche Breadth --##
    
    #- Calculate MVNH Breadth -#
      
      # Create csv to store results
      echo "total, tmax, ndvi, elev, cor, week, individual, scientificname, studyid, year" > ./out/niche_determinant_anthropause.csv

      # Make log file to track successful outputs
      echo "studyid, individual, scientificname, year, status, week" > out/niche_log.csv
  
      # Inputs: db:event_clean  + out filepath + no. cores
      # Outputs: csv 
      
      sbatch $src/hpc/run_calc_niche_breadth.sh
    #
  




  ##


####---- Inferential Models ----####

 #- Fit space use models -#
      
      # Inputs: space use csv + trait csv + no. cores
      # Outputs: model rdata objects 
       
      # SLURM

      sbatch $src/hpc/run_fit_space_use_dot_models.sh # dot
      sbatch $src/hpc/run_fit_space_use_interactive_models.sh # interactive
      sbatch $src/hpc/run_fit_space_use_additive_models.sh

      # ON DEMAND:
      # conda activate brms
      # Rscript $wd/analysis/src/workflow/fit-space-use-models.r $wd/out/dbbmm_size.csv $wd/out/single_species_models/area 24 10 10000 5

    #
  
    #- Generate model results and diagnostic sheets-#
      # TODO: Scott work scratch code into workflow
      # 
      # Only runs interactively:
      plot_space_use.r
    #

    #- Fit niche breadth models -#
      
      sbatch $src/hpc/run_fit_niche_breadth_dot_models.sh
      sbatch $src/hpc/run_fit_niche_breadth_additive_models.sh
      sbatch $src/hpc/run_fit_niche_breadth_interactive_models.sh #interaction model

  ##-- Intra-Individual Analysis --##
  
    #Area
    sbatch $src/hpc/run_fit_intra_ind_mod_additive_area.sh
    sbatch $src/hpc/run_fit_intra_ind_mod_int_area.sh
    
    #Niche
    sbatch $src/hpc/run_fit_intra_ind_mod_additive_niche.sh
    sbatch $src/hpc/run_fit_intra_ind_mod_int_niche.sh
  ##


  
####



####---- Outputs ----####

  ##-- Figures --##
  
  ##


  ##-- Other Outputs --##
  
  ##

####