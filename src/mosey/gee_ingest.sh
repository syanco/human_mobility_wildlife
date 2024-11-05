#!/bin/bash

# Overall Script Function:
# 1. evaluates the command run in wf-mosey_env.sh with args for where to find the db,
#     where to output the CSV's to, and where to send them within GEE
# 2. downloads trimmed event data from the processed_data dir in this repo to local CSVs, with some processing:
#       - renames event_id to anno_id
#       - formats the time stamp
#       - creates group numbers based on groupSize because of the limitation of GEE to authenticate a certain amount of data at once
# 3. if CSV has data, it is uploaded to GCS
# 4. uploads from GCS to GEE

#   Need to have gee environment activated

#TODO: could pass in optional parameters for upload gcs, import gee, delete csv
#TODO: still need to generalize to use any entity (population, study, etc.).
# Right now partially there but still coded to use population
#TODO: figure out how to activate pyenv gee inside the script
#TODO: figure out how to check if gee env is activated

eval "$(docopts -h - : "$@" <<EOF
Usage: gee_ingest.sh [options] <argv> ...
Options:
      --help     Show help options.
      --version  Print program version.
      --db=db    Path to the mosey database.
----
gee_ingest.sh 0.1
EOF
)"

#Fails w/in script
#setopt interactivecomments

#For testing:
# wd=~/projects/ms2/analysis/poc/mosey_env/mosey_env1
# mkdir -p $wd
# cd $wd
# sesnm=main
# geePtsP=users/benscarlson/projects/ms2/poc/mosey_env/mosey_env1
# gcsURL=gs://mol-playground/benc/projects/ms2/poc/mosey_env/mosey_env1/ingest
# db=~/projects/ms2/analysis/main/data/mosey.db
# outP=data/anno/ingest

# sesnm=${argv[0]}  # The name of the session
# geePtsP=${argv[1]} # The path to the gee folder that will hold the point assets
# gcsURL=${argv[2]} # The gcs url to the folder that will hold csvs for import to gee
# outP=${argv[3]} # The path to the folder that will hold csvs for import to gcs
sesnm=$1  # The name of the session
geePtsP=$2 # The path to the gee folder that will hold the point assets
gcsURL=$3 # The gcs url to the folder that will hold csvs for import to gee
outP=$4 # The path to the folder that will hold csvs for import to gcs

#Set defaults for optional paramters
[[ -z "$db" ]] && db=processed_data/mosey_mod.db
# [[ -z "$db" ]] && db=processed_data/mosey_mod.db

# Local parameters
groupSize=500000 #Pass in as optional argument
#table="forage_event" #TODO: probably need to pass in an sql statement or point to sql file? 
entity=study #Pass in as optional argument

#Get the session id
# sesid=$(sqlite3 $db "select ses_id from session where ses_name = '$sesnm' and table_name = 'study'")

mkdir -p $outP # should this be "gcsOutP" instead?
# earthengine create folder -p $geePtsP

# entIds=($(mlr --csv --opprint filter '$run == 1' then cut -f individual_id ctfs/$entity.csv | tail -n +2))
entIds=($(mlr --csv --opprint filter '$run == 1' then cut -f study_id ctfs/$entity.csv | tail -n +2))

n=${#entIds[@]}

echo Loading $n $entity.


# use for loop to read all values and indexes
for (( i=0; i<=${n}; i++ ));
do
  #entId=30
  #entName=east_portugal
  #i=1
  entId=${entIds[$i]}
  # entName=${names[$i]}
  
  echo "*******"
  echo "Start processing ($entity id: $entId)"
  echo "*******"
  
  csv=${outP}/${entId}.csv
  export csv=${outP}/${entId}.csv
  export gcsCSV=${gcsURL}/${entId}.csv
  export geePts=${geePtsP}/${entId}

  
  #Earthengine can't do annotation with tasks greater than 1e6 points
  #So need to break them up. row_number returns 1..n. Subtract 1 so that the
  # resulting number of groups will be correct. Since groupSize is an integer,
  # the result will be cast to an integer (equivilant to floor() operation)
  
  # Latest sql. Selects by pop instead of study
  # Calculates groups according to lon, lat, and also sorts the final dataset
  # although final sorting should not matter as long as group assignment is ordered
  sql="select f.event_id as anno_id, e.lon, e.lat,
    strftime('%Y-%m-%dT%H:%M:%SZ',e.timestamp) as timestamp, 
    	(row_number() over (order by e.lon, e.lat)-1)/$groupSize as grp 
    from event_trim f 
        inner join event e
    on f.event_id = e.event_id
    where f.study_id = ${entId}"
    
#     inner join event e on f.event_id = e.event_id
# 	  inner join study study_id on seg.pop_id = pop.pop_id
#     where pop.pop_id = $entId and pop.ses_id = $sesid
# 	order by anno_grp, lon, lat"

    # inner join forage_seg fs on f.fs_id = fs.fs_id 
    # inner join segment seg on fs.seg_id = seg.seg_id
 
  echo Extracting data...

  #echo $sql
  
  /usr/bin/env sqlite3 -header -csv $db "$sql;" > $csv
  
  #Check number of rows and skip if there are none
  #If no rows, file will be empty and rows=0
  #If there are rows, file will have header so num records is rows-1
  rows=$(  cat $csv | wc -l )
  
  if [ $rows = 0 ]; then
    echo "Extracted 0 rows, skipping."
    rm -f $csv
    continue
  fi
  
  echo "Extracted $(($rows-1)) rows" #subtract the header
  
  #---- Upload file to GCS
  echo Uploading file to gcs...
  gsutil -q cp -r $csv $gcsCSV

  #---- Import file into GEE
  echo Starting GEE import task...
  earthengine upload table  --asset_id $geePts $gcsCSV --x_column lon --y_column lat --force

  echo "Removing csv. $csv is: $csv"
  #---- Cleanup
  rm -f $csv
done

echo "Script complete"


#---- OLD CODE ----#

  #This is the old SQL I used to annotate the covid data and my schema
  # sql="select f.event_id as anno_id, f.study_id, e.lon, e.lat, 
  #     strftime('%Y-%m-%dT%H:%M:%SZ',e.timestamp) as timestamp,
  #     (row_number() over (order by f.study_id)-1)/${groupSize} anno_grp
  #   from ${table} f
  #   inner join event e
  #   on f.event_id = e.event_id
  #   where f.study_id = ${studyId}"

  #This is older sql
  #TODO: try ordering by lon,lat instead of event_id. This will tend to group
  # records geographically, which may increase performance. One task was
  # failing with memory use errors it seems b/c the wide extent of the points
  # Another approach would be to have GEE make geographic clusters.
  # sql="select f.event_id as anno_id, i.study_id, e.lon, e.lat, 
  #   strftime('%Y-%m-%dT%H:%M:%SZ',e.timestamp) as timestamp, 
  #   	(row_number() over (order by f.event_id)-1)/${groupSize} as anno_grp 
  #   from forage_event f 
  #   inner join event e on f.event_id = e.event_id
  #   inner join forage_seg fs on f.fs_id = fs.fs_id 
  #   inner join segment seg on fs.seg_id = seg.seg_id
  #   inner join individual i on seg.individual_id = i.individual_id
  #   where i.study_id = ${studyId} and fs.ses_id = ${sesid}"
  
#---- Timing various sql statements ----#

#SQL 1
  # sql="select f.event_id, f.study_id, e.lon, e.lat, strftime('%Y-%m-%dT%H:%M:%SZ',e.timestamp) as timestamp \
  #   from event_forage f \
  #   inner join ( \
  #     select event_id, lon, lat, timestamp from event where study_id = ${studyId} \
  #   ) e on f.event_id = e.event_id \
  #   where f.study_id = ${studyId}"

#Note tests with 3 million event dataset show this method is just as fast!
#Test timing when writing the full dataset to csv

#SQL 2
  # sql="select f.event_id, f.study_id, e.lon, e.lat, strftime('%Y-%m-%dT%H:%M:%SZ',e.timestamp) as timestamp
  #   from event_forage f
  #   inner join event e
  #   on f.event_id = e.event_id
  #   where f.study_id = ${studyId}"

#SQL 3
  # sql="select f.event_id, f.study_id, e.lon, e.lat, strftime('%Y-%m-%dT%H:%M:%SZ',e.timestamp) as timestamp
  #   from event_forage f
  #   inner join event e
  #   on f.event_id = e.event_id
  #   where f.study_id = ${studyId} and e.study_id = ${studyId}"
  
#SQL 2 is simplist and is just as fast so use that query.
  #10763606 SQL 1 = 0.75 real
  #10763606 SQL 2 = 0.75 real
  #8863543 SQL 1 = 18.69 real
  #8863543 SQL 2 = 8.90 real
  #8863543 SQL 1 = 9.40 real
  #8863543 SQL 3 = 9.18 real
