# Master workflow script for annotating database event table with environmental
# layers via Google Earth Engine (GEE).
# Notes:
# - Before executing, this script requries 
#   the user to create a GEE account and bucket, create subdirectories,
#   change filepaths, and create a local python environment.
# - Using the rgee package with a local python environment requires user
#   to configure a python env path in the rgee config file
# - Each shell script that is launched from this script needs to complete
#   running prior to launching the next.

wd=/Users/juliet/Documents/OliverLab/covid_paper/repositories/human_mobility_wildlife
src=$wd/src
db=$wd/processed_data/mosey_mod.db
export MOSEYENV_SRC=/Users/juliet/Documents/OliverLab/covid_paper/repositories/human_mobility_wildlife/src/mosey
export RETICULATE_PYTHON=/opt/anaconda3/envs/hmw_py3-9/bin/python
cd $wd

sql="select study_id, 1 as run
from study_trim
order by study_id"

sqlite3 -header -csv $db "$sql;" > ctfs/study.csv

# name of GEE project to hold 
export geePtsP=projects/covid-mvmnt-2024-440720/assets/tracks
# local folder that holds the csv files to be ingested into gee
export csvP=/Users/juliet/Documents/OliverLab/covid_paper/repositories/human_mobility_wildlife/gee_data/csvs_gee_ingest
# communicates which environmental data to annotate with in GEE
export envP=/Users/juliet/Documents/OliverLab/covid_paper/repositories/human_mobility_wildlife/ctfs/env.csv
# local folder that holds the annotated CSV files after GEE step is complete
export annoP=/Users/juliet/Documents/OliverLab/covid_paper/repositories/human_mobility_wildlife/gee_data/annotated

##= -- GCS & GEE DIRS --
export gcsBucket=covid-mvmnt-2024
# dir for CSVs in GCS that will be imported GEE
export gcsInP=ingest_gee
# output folder for annotated CSVs (excluding bucket name) in GCS
export gcsOutP=annotated

# -- GEE URLs --
export gcsInURL=gs://${gcsBucket}/${gcsInP} #This is the url to the gee ingest folder
export gcsOutURL=gs://${gcsBucket}/${gcsOutP} #This is the url to the output folder (includes bucket)

# send csv files to GCS & GEE ingest_GEE dir (launches 1 job per species dataset)
$MOSEYENV_SRC/gee_ingest.sh trial_1 $geePtsP $gcsInURL $csvP

# generate annotations in GEE (launches 3 jobs per species dataset)
$MOSEYENV_SRC/mosey_anno_gee.sh $geePtsP $gcsOutP 

# import annotated data into mosey DB:
# first add table and columns to db to receive annos for 3 env layers
# (simply creating the database infrastructure to be populated)
sqlite3 $db "alter table event_trim add column tmax REAL;"
sqlite3 $db "alter table event_trim add column ndvi REAL;"
sqlite3 $db "alter table event_trim add column elev REAL;"
# populate the database with the data:
# points and annotations are stored in event_clean
# for testing, see db/anno_test.sql for ddl to create event_test
$MOSEYENV_SRC/import_anno.sh $gcsOutURL $annoP $db --table event_trim
