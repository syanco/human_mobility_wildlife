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
#TODO: figure out how to activate env inside the script
#TODO: figure out how to check if gee env is activated without manual step

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

# args:
# sesnm=${argv[0]}  # The name of the session
# geePtsP=${argv[1]} # The path to the gee folder that will hold the point assets
# gcsURL=${argv[2]} # The gcs url to the folder that will hold csvs for import to gee
# outP=${argv[3]} # The path to the folder that will hold csvs for import to gcs

# use bash syntax instead of zsh
sesnm=$1  # The name of the session
geePtsP=$2 # The path to the gee folder that will hold the point assets
gcsURL=$3 # The gcs url to the folder that will hold csvs for import to gee
outP=$4 # The path to the folder that will hold csvs for import to gcs

# Set default
[[ -z "$db" ]] && db=processed_data/mosey_mod.db

# Local parameters
groupSize=500000 #Pass in as optional argument
#table="forage_event" #TODO: probably need to pass in an sql statement or point to sql file? 
entity=study #Pass in as optional argument

# Get the session id
# sesid=$(sqlite3 $db "select ses_id from session where ses_name = '$sesnm' and table_name = 'study'")

mkdir -p $outP
# earthengine create folder -p $geePtsP

entIds=($(mlr --csv --opprint filter '$run == 1' then cut -f study_id ctfs/$entity.csv | tail -n +2))

echo "Extracted study IDs:"
printf '%s\n' "${entIds[@]}"

n=${#entIds[@]}

echo Loading $n $entity.


# use for loop to read all values and indexes
for (( i=0; i<${n}; i++ ));
do

  entId=${entIds[$i]}
  
  echo "*******"
  echo "Start processing ($entity id: $entId)"
  echo "*******"
  
  csv=${outP}/${entId}.csv
  export csv=${outP}/${entId}.csv
  export gcsCSV=${gcsURL}/${entId}.csv
  export geePts=${geePtsP}/${entId}

  
  # Earthengine can't do annotation with tasks greater than 1e6 points
  # So need to break them up. row_number returns 1..n. Subtract 1 so that the
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
 
  echo Extracting data...
  
  /usr/bin/env sqlite3 -header -csv $db "$sql;" > $csv
  
  # Check number of rows and skip if there are none
  # If no rows, file will be empty and rows=0
  # If there are rows, file will have header so num records is rows-1
  rows=$(  cat $csv | wc -l )
  
  if [ $rows = 0 ]; then
    echo "Extracted 0 rows, skipping."
    rm -f $csv
    continue
  fi
  
  echo "Extracted $(($rows-1)) rows" # subtract the header
  
  #---- Upload file to GCS
  echo Uploading file to gcs...
  gsutil -q cp -r $csv $gcsCSV

  #---- Import file into GEE
  echo Starting GEE import task...
  earthengine upload table  --asset_id $geePts $gcsCSV --x_column lon --y_column lat --force

  #---- Cleanup
  rm -f $csv
done

echo "Script complete; GEE ingest job has launched."
