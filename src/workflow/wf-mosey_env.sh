# Annotate movement data using mosey env.
# 
# MUST BE RUN INTERACTIVELY 
# (b/c progress depends on remote processing times on google.)
#
#
# #---- Get the database from hpc ----#
# 
# #Connect to vpn
# /opt/cisco/anyconnect/bin/vpn connect access.yale.edu
# # 
# wd=/Users/scottyanco/Documents/covid-19_movement/out/anno
# # dbr=/home/sy522/project/covid-19_movement/processed_data/mosey_swap_mod.db
dbr=/scratch/julietcohen/covid_movement/human_mobility_wildlife/processed_data/mosey_mod.db
# 
# cd $wd

# #check file size
# ssh grace "ls -lh $dbr" #8.3GB
# 
# #download the file
# scp grace:$dbr mosey_mod.db
# 
# #disconnect from vpn
# /opt/cisco/anyconnect/bin/vpn disconnect

#----
#---- mosey_env
#----


wd=/Users/juliet/Documents/OliverLab/covid_paper/repositories/human_mobility_wildlife
src=$wd/src
# db=$wd/processed_data/mosey_mod.db
db=$wd/processed_data/mosey_mod.db

# define filepath to lcoal repo's GEE scripts
export MOSEYENV_SRC=/Users/juliet/Documents/OliverLab/covid_paper/repositories/human_mobility_wildlife/src/mosey #for mosey_anno_gee.sh

cd $wd

#--
#-- Create study.csv file from the trimmed table in processed_data/mosey_mod.db
#--

# mkdir ctfs

sql="select study_id, 1 as run
from study_trim
order by study_id"

sqlite3 -header -csv $db "$sql;" > ctfs/study.csv

#--
#-- Set up variables
#--

geePtsP=/Users/juliet/Documents/OliverLab/covid_paper/repositories/human_mobility_wildlife/processed_data #folder holding the gee point datasets
# geePtsP=project/covid-mvmnt/assets/tracks #folder holding the gee point datasets
gcsBucket=covid-mvmnt-2024
gcsInP=ingest_gee #This holds the CSVs that will be imported into GEE
gcsOutP=annotated #This is the output folder for annotated csvs (excluding bucket) in GCS
csvP=/Users/juliet/Documents/OliverLab/covid_paper/repositories/human_mobility_wildlife/csvs_for_gee_ingest #local folder that holds the csv files to be ingested into gee
annoP=data/anno/annotated #local folder that holds the annotated csv files
envP=analysis/ctfs/env.csv

#TODO: don't require sesid (session ID)
#sesid=full_wf

gcsInURL=gs://${gcsBucket}/${gcsInP} #This is the url to the gee ingest folder
gcsOutURL=gs://${gcsBucket}/${gcsOutP} #This is the url to the output folder (includes bucket)

#NASA/ORNL/DAYMET_V4; Tmax: band 4. (daymet_tmax). Daily with start/end date
#NASA/ORNL/DAYMET_V4; Tmin: band 5. (daymet_tmin). Daily with start/end date
#MODIS/006/MOD11A1; LST_day_1km: band 0. (lst). Daily with start/end date
#CSP/HM/GlobalHumanModification; gHM: band 0. (ghm). GHM is a strange layer. A single, static image in an image collection. No start/end dates.
# projects/map-of-life/diegoes/dist2road_USA_full; distance: band 0. (dist2road)
#MODIS/MOD09GA_006_NDVI; NDVI: band 0. (ndvi). Daily with start/end date. 
#CGIAR/SRTM90_V4; elevation: band 0. (elev). Static image. Has start/end date but this is for a 10 day window in Feb 2000

# Studies for testing
# study_id, num
# 1450604622,	3067
# 619097045,	4469127

#----
#---- Update the database structure ----#
#----

# See db/edit_structure.sql

#----
#---- Import studies locally, modify, then send to bucket & GEE
#----

## analysis/ - let ths be passed by arg
$MOSEYENV_SRC/gee_ingest.sh trial_1 $geePtsP $gcsInURL $csvP #$sesid

# #----
# #---- Annotate 
# #----
# 
# #Don't run this until all import tasks have finished
# #https://code.earthengine.google.com/tasks
# 
# $MOSEYENV_SRC/mosey_anno_gee.sh $geePtsP $gcsOutP #$envP #"${envs[*]}"
# 
# #----
# #---- Import into mosey ----#
# #----
# 
# # add columns to db to receive annos
# sqlite3 $db "alter table event_trim add column tmax REAL;"
# # sqlite3 $db "alter table event_clean add column tmin REAL;"
# # sqlite3 $db "alter table event_clean add column lst REAL;"
# sqlite3 $db "alter table event_trim add column ndvi REAL;"
# sqlite3 $db "alter table event_trim add column elev REAL;"
# # sqlite3 $db "alter table event_clean add column dist2road REAL;"
# 
# #points and annotations are stored in event_clean
# #for testing, see db/anno_test.sql for ddl to create event_test
# $MOSEYENV_SRC/import_anno.sh $gcsOutURL $annoP $db --table event_trim
# 
# # Send db back to HPC
# scp $db grace:$dbr
