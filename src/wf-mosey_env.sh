#---- Get the database from hpc ----#

#Connect to vpn
/opt/cisco/anyconnect/bin/vpn connect access.yale.edu

wd=analysis/anno
dbr=/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod.db

cd $wd

#check file size
ssh grace "ls -lh $dbr" #8.3GB

#download the file
scp grace:$dbr data

#disconnect from vpn
/opt/cisco/anyconnect/bin/vpn disconnect

#----
#---- mosey_env
#----

pyenv activate gee

wd=~/projects/covid/analysis/anno
src=~/projects/covid/src
db=data/mosey_mod.db

export MOSEYENV_SRC=~/projects/mosey_env/src #for mosey_anno_gee.sh

cd $wd

#--
#-- Create study.csv file since I don't have one
#--

mkdir ctfs

sql="select study_id, study_name, 1 as run
from study
where study_id in (select distinct study_id from event_clean)
order by study_id"

/usr/bin/time sqlite3 -header -csv $db "$sql;" > ctfs/study.csv

#--
#-- Set up variables
#--

geePtsP=users/benscarlson/projects/covid/tracks #folder holding the gee point datasets
gcsBucket=mol-playground
gcsInP=benc/projects/covid/ingest_gee #This holds the csvs that will be imported to gee
gcsOutP=benc/projects/covid/annotated #This is the output folder for annotated csvs (excluding bucket)
csvP=data/anno/ingest_gee #local folder that holds the csv files to be ingested into gee
annoP=data/anno/annotated #local folder that holds the annotated csv files

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

#TODO: db should be optional
/usr/bin/time $MOSEYENV_SRC/main/gee_event_import.sh $geePtsP $gcsInURL $db $csvP #$sesid

#----
#---- Annotate 
#----

#Don't run this until all import tasks have finished
#https://code.earthengine.google.com/tasks

/usr/bin/time $MOSEYENV_SRC/main/mosey_anno_gee.sh $geePtsP $gcsOutP #"${envs[*]}"

#----
#---- Import into mosey ----#
#----

#points and annotations are stored in event_clean
#for testing, see db/anno_test.sql for ddl to create event_test
$MOSEYENV_SRC/main/import_anno.sh $gcsOutURL $annoP $db --table event_clean

#----
#---- Upload the database----#
#----

#Connect to vpn
/opt/cisco/anyconnect/bin/vpn connect access.yale.edu

wd=analysis/anno
dbAnno=/gpfs/loomis/project/jetz/sy522/covid-19_movement/processed_data/mosey_mod_anno.db

cd $wd

#upload the file (also rename it)
scp data/mosey_mod.db grace:$dbAnno 

#disconnect from vpn
/opt/cisco/anyconnect/bin/vpn disconnect