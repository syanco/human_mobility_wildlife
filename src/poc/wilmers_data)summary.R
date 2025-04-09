library(lubridate)
stdid <- 1891257164
x <- evt_trm %>% 
  filter(study_id == stdid) %>% 
  group_by(individual_id) %>% 
  mutate(ts = ymd_hms(timestamp),
    dur = ts-dplyr::lag(ts)) %>% 
  summarise(n = n(),
            interval = median(dur, na.rm = T),
            interval_min = interval/60,
            interval_hrs = interval_min/60)
x
