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

# Top-level workfow script to execute/describe steps to run analyses associated 
# with the North American COVID-19 & Animal Movemnt project. 

# This project is led by: Ruth Oliver, Diego Ellis-Soto, Scott Yanco, Brett 
# Jesmer, & Walter Jetz (PI). 

# This workflow is generally built to be run in a high performance computing 
# (HPC) environment. Thus, most steps in the workflow simply use `sbatch` to 
# submit a job to the slurm manager.  However, we also include (commented) code 
# to run each step without submission to a slurm manager if desired 
# (not recommended).

# This workflow also assumes a specific directory strcuture.  See README for 
# detail son directory strcuture.

# CONDA ENVIRONMENTS
# 
# Much of this code is meant to be executed within coda environments. The 
# directory /analysis/src/conda_envs contains .yml files from which the 
# requisite conda environments can be built. See readme for details on building 
# conda environments

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#

  # TODO:
  
  #   * relativize all sbatch .sh scripts (i.e., call .r files relative to $wd)
  #   * Fix breezy paths (e.g., i think the `rd(src/startup.r)` call won't work)

#------------------------------------------------------------------------------#

####----  Initialization  ----####
  
  # define working directory
  wd=~/project/covid-19_movement
  src=$wd/analysis/src
  
  # Go to working directory
  cd $wd
  
  # Make scripts executable (mostly to make docopts visible via `-h` flag)
  chmod +x $src/workflow/clean_movement.r
  chmod +x $src/workflow/compute-cbg-area.r
  chmod +x $src/workflow/annotate-events-cbg.r
  
  # Turn on miniconda (only necessary on HPC when not using SLURM)
  # module load miniconda

####



####----  Prep and Clean Data ----####

  ##-- Working DB --##

    # Make a copy of db to be modified
    cp $wd/raw_data/mosey.db $wd/processed_data/mosey_mod.db
  ##

  ##-- Run Cleaning --##
  
    # Inputs:   db:event + analysis/ctfs/dates.csv
    # Outputs:  db:event_mod + db:event_clean
    # Depends:  clone of mosey_env (https://github.com/benscarlson/mosey_env)
  
    # SLURM:
    sbatch $src/hpc/run_clean_movement.sh

    # ON DEMAND:
    # conda activate spatial
    # Rscript $src/workflow/clean_movement.r --db $wd/processed_data/mosey_mod.db
  
    # TODO: 
    #   * Remove hardcoded study filters
    #   * Simplify date filtering
    #   * Add spatial filtering (events outside US).
    #   * Allow quantile thresholds to be passed by arg
  ##

####



####----  Annotate Data ----####

  ##-- Environmental Annotations --#
  
    # Inputs: db:event_clean + analysis/ctfs/env.csv
    # Outputs: db:event_clean (annotated)
    
    # INTERACTIVE (script must be run interactively)
    $src/wf-mosey_env.sh
    
    # TODO:
    #   * Add conda environment for mosey_env
  
  ##
    
    
  ##-- Census Annotations --##
  
    #- Intersect events with cbg geometries -#
    
      # Inputs:   db:event_clean + cbg shp file
      # Outputs:  csv per job (event_id + cbg info)

      # DSQ: (Yale-specific)
      module load R/4.1.0-foss-2020b
      Rscript $src/workflow/create_intersection_joblist.r

      module load dSQ
      dsq --job-file $src/workflow/joblist.txt --mem-per-cpu 40g -t 2-

      # UPDATE WITH DATE
      sbatch dsq-joblist-2022-05-12.sh
      
      # TODO: 
      #   * Missing code for last step above ("dsq-joblist-2022-05-12.sh")
      #   * Rewrite for non-DSQ (for reproducibility)?
      #   * make conda env for this step, remove module r, dsq?
    
    #
    
    #- Compute cbg area -#

      # Inputs: cbg shp files
      # Outputs: csv (cbg info + area)

      # SLURM:
      sbatch $src/hpc/run_compute_cbg_area.sh
      
      # ON DEMAND:
      # Rscript $src/workflow/compute-cbg-area.r
      
      # TODO:
      #   * Fix/update docopts in compute-cbg-area.r
      #   * Pass .shp location as arg?
      #   * Use conda rather than module r?
    
    #
    
    #- Annotate events with cbg info -#
    
      # Inputs: db:event_clean + cbg intersection csv + cbg area csv
      # Outputs: csv (event_id + cbg info + cbg area)

      # SLURM:
      sbatch $src/hpc/run_annotate_events_cbg.sh
      
      # ON DEMAND:
      # Rscript $src/workflow/annotate-events-cbg.r

      # TODO:
      #   * Fix/update docopts in annotate-events-cbg.r
      #   * Pass input/output paths as arg?
      #   * Use conda rather than module r?
    
    #
  
  ##
  
  
  ##-- SafeGraph Annotations --##
  
    #- Process SafeGraph data -#
      
      # Inputs: safegraph txt files (one file per cbg/week)
      # Outputs: sg data csv (one file per cbg/week)

      # SLURM:
      sbatch $src/hpc/run_process_safegraph_data.sh
      
      # ON DEMAND:
      # Rscript $src/workflow/process-safegraph-data.r
      
      # TODO:
      #   * Fix/update docopts in process-safegraph-data.r
      #   * Pass input/output paths as arg?
      #   * Use conda rather than module r?
      #   * comment out/delete deprecated code line 118 and beyond?
    
    #

    #- Annotate events with SafeGraph -#
    
      # Inputs: db:event_clean + cbg info csv + sg data csv
      # Outputs: csv (event_id + timestamp + cbg info = cbg area + sg count)
  
      # SLURM:
      sbatch $src/hpc/run_annotate_events_safegraph.sh
      
      # ON DEMAND:
      # Rscript $src/workflow/annotate-events-safegraph.r
      
      # TODO:
      #   * Fix/update docopts in annotate-events-safegraph.r
      #   * Pass input/output paths as arg?
      #   * Use conda rather than module r?
  
    #

  ##
  
  
  ##-- Human Modification Annotations --#
  
    #- Annotate events with TNC GHM -#
    
      # Inputs: db:event_clean + ghm raster
      # Ouputs: csv (event_id + ghm)

      # SLURM:
      sbatch $src/hpc/run_annotate_events_ghm.sh
    
      # ON DEMAND:
      # Rscript $src/workflow/annotate-events-ghm.r
      
      # TODO:
      #   * Fix/update docopts in annotate-events-ghm.r
      #   * Pass input/output paths as arg?
      #   * Use conda rather than module r?
    
    #
  
  ##

####


####----  Analysis ----####

  ###---  Fit dBBMMs  ---###

# Activate covid env
# conda activate covid

# Make dir to hold fitted dBBMMs (only run once)
mkdir $wd/out/dbbmms

# Make log file to track successful outputs (only run once)
echo "species, ind_id, study_id, year, out_type, filename, produced, out_date" > out/dbbmm_log.csv

# Make big mem log file to track ind-year combos saved for the big mem parition
echo "species, ind_id, study_id, year, out_type, filename, produced, out_date" > out/dbbmm_bigmem_log.csv

# CHOOSE ONE:
  # Fit dBBMMs (run on local - I don not recommend)
  # Rscript $wd/analysis/src/0X-fit_dBBMMs.r -nc 2
  
  # Fit dBBMMs (run on local - I don not recommend)
  sbatch ~/project/covid-19_movement/analysis/src/hpc/submit_dBBMM.sh
  
# Calc dbbmm areas
echo "species, ind_id, study_id, year, wk, area, sg, ghm, cbg_area, pop, ndvi, lst, n, a_bb, fixmed, m_error" > ./out/dbbmm_size.csv

# CHOOSE ONE:
  # Calc sizes on HPC 
  sbatch ~/project/covid-19_movement/analysis/src/hpc/submit_size.sh
