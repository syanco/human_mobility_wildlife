# Annotate movement data using mosey env.
# 
# MUST BE RUN INTERACTIVELY 
# (b/c progress depends on remote processing times on google.)
#
#

#----
#---- mosey_env
#----

# pyenv activate gee
mamba activate annotate

wd=~/projects/covid-19_movement
src=$wd/analysis/src
db=processed_data/mosey_mod_anno.db

export MOSEYENV_SRC=~/projects/mosey_env/main #for mosey_anno_gee.sh

cd $wd

# #--
# #-- Create study.csv file since I don't have one
# #--
# 
# mkdir ctfs
# 
# sql="select study_id, study_name, 1 as run
# from study
# where study_id in (select distinct study_id from event_clean)
# order by study_id"
# 
# /usr/bin/time sqlite3 -header -csv $db "$sql;" > ctfs/study.csv

#--
#-- Set up variables
#--

geePtsP=users/syanco/covid-mvmnt/bg-pts #folder holding the gee point datasets
gcsBucket=covid-mvmnt-bucket
gcsInP=ingest_gee #This holds the csvs that will be imported to gee
gcsOutP=bg_annotated #This is the output folder for annotated csvs (excluding bucket)
csvP=out/ssf-background-pts/individual-files #local folder that holds the csv files to be ingested into gee
annoP=out/ssf-background-pts/annotated #local folder that holds the annotated csv files
envP=analysis/ctfs/env.csv

gcsInURL=gs://${gcsBucket}/${gcsInP} #This is the url to the gee ingest folder
gcsOutURL=gs://${gcsBucket}/${gcsOutP} #This is the url to the output folder (includes bucket)


#----
#---- Import studies into GEE
#----

# Load vector of inds in from log file
# Use miller to filter by run column and then take the study_id field
# need to use tail to remove first line, which is the header
# study ids have \r suffix. Need to remove these
indIds=($(mlr --csv --opprint cut -f ind out/ssf-background-pts/bg-log.csv | tail -n +2))
indIds=${indIds[@]%$'\r'} # remove \r suffix

# cd $csvP

# Loop thorugh inds to upload to gcs and gee
for IND in $indIds; do 

  gcsCSV=$gcsInURL/${IND}.csv
  geePts=$geePtsP/$IND

  #---- Upload file to GCS
  echo Uploading $IND to gcs...
  gsutil -q cp -r $csvP/${IND}.csv $gcsCSV

  #---- Import file into GEE
  echo Starting GEE import task for individual $IND...
  
  earthengine upload table $gcsCSV --asset_id $geePts --x_column lon --y_column lat --force
done


#----
#---- Annotate 
#----

#Don't run this until all import tasks have finished
#https://code.earthengine.google.com/tasks

envs=($(mlr --csv --opprint filter '$run == 1' then cut -f "env_id" $envP | tail -n +2))
bands=($(mlr --csv --opprint filter '$run == 1' then cut -f band $envP | tail -n +2))
colnames=($(mlr --csv --opprint filter '$run == 1' then cut -f col_name $envP | tail -n +2))

# Remove \r suffix
envs=(${envs[@]%$'\r'})
bands=(${bands[@]%$'\r'})
colnames=(${colnames[@]%$'\r'})

# echo ${envs[@]}
# # envs=${envs[@]:1}
# echo Annotating ${#indIds[@]} background sets.
# 
# # Find index of last individual sent to GEE, then subset the vector after that 
# # to restart loop if it stops before completion
# # indIDs=(${indIds[0]})
# value=1153545169
# 
# for i in "${!indIds[@]}"; do
#    if [[ "${indIds[$i]}" = "${value}" ]]; then
#        echo "${i}";
#    fi
# done
# 
# indIds=${indIds[@]:699}

for indId in $indIds; do 
  echo "Start processing individual ${indId}"

  points=$geePtsP/$indId
  
  # get length of an array
  n=${#envs[@]}

  # use for loop to read all values and indexes
  for (( i=0; i<${n}; i++ ));
  do
  
    #i=0
    
    #TODO: do this as default if user doesn't pass in col_name info
    #envN=${env##*/} #gets the name (w/o path) of the env variable
    
    #TODO: check to see if $points exists in gee before annotating
    # earthengine asset info $points
    # earthengine asset info x
    # return_value=$?
    
    #TODO: need to handle band, colname as optional parameters
    # if column is not present don't pass parameters

    #echo "index: $i, env: ${envs[$i]}, band: ${bands[$i]}, col name: ${colnames[$i]}"
    out=$gcsOutP/${indId}_${colnames[$i]} #do not include url, bucket, or file extension
    
    echo Annotating "env: ${envs[$i]}, band: ${bands[$i]}, col name: ${colnames[$i]}"
    $MOSEYENV_SRC/anno_gee.r $points ${envs[$i]} $out -b ${bands[$i]} -c ${colnames[$i]}
    
  done

done

#----
#---- Download annotations from google ----#
#----

for IND in ${indIds}; do
	echo "Downloading individual ${IND}"

  # get length of an array
  n=${#envs[@]}

  # use for loop to read all values and indexes
  for (( i=0; i<${n}; i++ )); 
  do
  
  	# for env in "${envs[@]}"
	  # do
  	#  envN=${env##*/} #gets the name (w/o path) of the env variable
	  #  echo "Importing $envN"
	  #   annoN=${studyId}_${envN}
		
		
    #i=0
    
		#env=users/benscarlson/projects/ms3/dist2forest

    echo "Importing ${colnames[$i]}"
    
		annoN=${IND}_${colnames[$i]}
		annoPF=$annoP/${annoN}.csv
		gcsCSV=$gcsOutURL/${annoN}.csv

		#check here if file exists and skip if it doesn't exist
		#https://stackoverflow.com/questions/48676712/how-to-check-if-any-given-object-exist-in-google-cloud-storage-bucket-through-ba
		
		#gsutil ls $gcsOutURL/${annoN}_*.csv
		gsutil -q stat $gcsOutURL/${annoN}_*.csv

		return_value=$? #returns 0 if files exist, 1 if there are no results

		if [ $return_value = 1 ]; then
			echo "$gcsCSV does not exist. Skipping."
			continue
		fi
		
		echo "Downloading $gcsOutURL/${annoN}_*.csv ..."
		
		gsutil cp $gcsOutURL/${annoN}_*.csv $annoP
	done
done
		
