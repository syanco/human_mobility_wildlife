cd raw_data/
scp mosey.db grace:~/project/covid-19_movement/raw_data/
scp mosey_swap.db grace:~/project/covid-19_movement/raw_data/
cp mosey_mod.db mosey_merged.db  
scp mosey_merged.db grace:~/project/covid-19_movement/raw_data/
  
cd ../processed_data
scp mosey_mod.db grace:~/project/covid-19_movement/raw_data/
scp mosey_swap.db grace:~/project/covid-19_movement/raw_data/