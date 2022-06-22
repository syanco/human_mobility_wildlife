evt <- tbl(db,'event_clean')
indtb <- tbl(db, "individual") 
stdtb <- tbl(db, 'study')


inds <- evt %>% 
  pull(individual_id) %>% 
  # collect() %>% 
  unique()

studies <- indtb %>% 
  filter(individual_id %in% inds) %>% 
  pull(study_id) %>% 
  unique()

foc_std <- stdtb %>% 
  filter(study_id %in% studies) %>% 
  collect()

write_csv(foc_std, "./out/study_table.csv")
