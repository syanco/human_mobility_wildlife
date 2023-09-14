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
# 
wd=/Users/scottyanco/Documents/covid-19_movement/out/anno
# dbr=/home/sy522/project/covid-19_movement/processed_data/mosey_swap_mod.db
dbr=/home/sy522/project/covid-19_movement/processed_datamosey_mod_2023.db

cd $wd

#check file size
ssh grace "ls -lh $dbr" #8.3GB

#download the file
scp grace:$dbr mosey_mod.db
# 
# #disconnect from vpn
# /opt/cisco/anyconnect/bin/vpn disconnect

#----
#---- mosey_env
#----


wd=/Users/scottyanco/Documents/covid-19_movement
src=$wd/analysis/src
# db=$wd/processed_data/mosey_mod.db
db=$wd/processed_data/mosey_mod_2023.db

export MOSEYENV_SRC=/Users/scottyanco/Documents/covid-19_movement/analysis/src/mosey #for mosey_anno_gee.sh

cd $wd

#--
#-- Create study.csv file since I don't have one
#--

# mkdir ctfs

sql="select study_id, 1 as run
from study_trim
order by study_id"

sqlite3 -header -csv $db "$sql;" > analysis/ctfs/study.csv

#--
#-- Set up variables
#--

geePtsP=users/syanco/covid-mvmnt/tracks_2023 #folder holding the gee point datasets
# geePtsP=project/covid-mvmnt/assets/tracks #folder holding the gee point datasets
gcsBucket=covid-mvmnt-bucket
gcsInP=ingest_gee #This holds the csvs that will be imported to gee
gcsOutP=annotated #This is the output folder for annotated csvs (excluding bucket)
csvP=out/anno/individual-files #local folder that holds the csv files to be ingested into gee
annoP=data/anno/annotated #local folder that holds the annotated csv files
envP=analysis/ctfs/env.csv

#TODO: don't require sesid
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
#---- Import studies into GEE
#----


## analysis/ - let ths be passed by arg
$MOSEYENV_SRC/gee_ingest.sh re_anno $geePtsP $gcsInURL $csvP #$sesid

#----
#---- Annotate 
#----

#Don't run this until all import tasks have finished
#https://code.earthengine.google.com/tasks

$MOSEYENV_SRC/mosey_anno_gee.sh $geePtsP $gcsOutP #$envP #"${envs[*]}"

#----
#---- Import into mosey ----#
#----

# add columns to db to receive annos
sqlite3 $db "alter table event_trim add column tmax REAL;"
# sqlite3 $db "alter table event_clean add column tmin REAL;"
# sqlite3 $db "alter table event_clean add column lst REAL;"
sqlite3 $db "alter table event_trim add column ndvi REAL;"
sqlite3 $db "alter table event_trim add column elev REAL;"
# sqlite3 $db "alter table event_clean add column dist2road REAL;"

#points and annotations are stored in event_clean
#for testing, see db/anno_test.sql for ddl to create event_test
$MOSEYENV_SRC/import_anno.sh $gcsOutURL $annoP $db --table event_trim

# Send db back to HPC
scp $db grace:$dbr
