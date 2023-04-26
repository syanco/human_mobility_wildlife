##-- Step Selection Function --##

#- Generate background points
# Inputs: db:event_clean
# Outputs: individual csvs 

# create log file which is then used as ctf for annotation step
echo "ind, run" > ./out/ssf-background-pts/bg-log.csv

# SLURM:
sbatch $src/hpc/run_generate_background_points.sh

# ON DEMAND:
# Rscript $src/workflow/generate-background-points.r
#


#- Environmental annotations for background points 

# Inputs: Background csvs + analysis/ctfs/env.csv
# Outputs: csv per individual X variable

# INTERACTIVE - RUN LOCAL (script must be run interactively and not on HPC...)

#pull bg files to local
scp -r grace:~/project/covid-19_movement/out/ssf-background-pts ~/projects/covid-19_movement/out
$src/workflow/wf-mosey_env-BG.sh
scp -r out/ssf-background-pts/annotated grace:~/project/covid-19_movement/out/ssf-background-pts/annotated
# TODO:
#   * Add conda environment for mosey_env

##

#- SafeGraph and GHM annotations for background pts

# Inputs: Background csvs
# Outputs: csv per job (event_id + cbg info)

# # module load R/4.1.0-foss-2020b
# conda activate covid
# # Rscript $src/workflow/create_bg_annotation_joblist.r
# Rscript $src/workflow/create_bg_annotation_joblist_moose.r
# 
# module load dSQ
# dsq --job-file $src/hpc/annotation-joblist.txt --mem-per-cpu 100g -t 2:00:00 -p pi_jetz
# 
# # UPDATE WITH DATE
# sbatch dsq-annotation-joblist-2022-08-10.sh

# ALTERNATIVE WAY USING FOREACH
# The above method keeps runign out of memore - not sure exactly why
# I re-wrote the scrip using MC parallelization and a submit script
# that calls on big mem so we cna have 100GB per core
sbatch $src/hpc/run_annotate_background_ghm_sg_ALL.sh

##

#- Fit SSFs
# Inputs: 
# Outputs:  

# SLURM:
sbatch $src/hpc/run_fit_SSF_models.sh

# ON DEMAND:
# Rscript $src/workflow/generate-background-points.r
#