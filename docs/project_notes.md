## TODO  
* Export FINAL conda env .yml files  
  * covid  
  * spatial  
* Clean comments on script 01  
* Check that dBBMM log file count matches the out/ directory


## Notes  

*Fitting dBBMMs for each ind x yr combination.

*I'm getting an error in thr prep and clean script on the HPC - that I DO NOT GET when run local:  
```
Error: Problem with `mutate()` column `sl`.
ℹ `sl = distGeo(cbind(lon, lat), cbind(lag_lon, lag_lat))`.
✖ longitude > 360
ℹ The error occurred in group 1068: individual_id = 1891176543.
Backtrace:
     █
  1. ├─`%>%`(...)
  2. ├─dplyr::mutate(...)
  3. ├─dplyr:::mutate.data.frame(...)
  4. │ └─dplyr:::mutate_cols(.data, ..., caller_env = caller_env())
  5. │   ├─base::withCallingHandlers(...)
  6. │   └─mask$eval_all_mutate(quo)
  7. ├─geosphere::distGeo(cbind(lon, lat), cbind(lag_lon, lag_lat))
  8. │ └─geosphere:::.pointsToMatrix(p1)
  9. │   └─base::stop("longitude > 360")
 10. └─base::.handleSimpleError(...)
 11.   └─dplyr:::h(simpleError(msg, call))
Execution halted
```
For now I'm just scp'ing the processed db to the HPC to run the dbbmms

* For the HR size - consider just a lmm of size ~ trt X yr + (1|ind) as basic model

### General

* Proposed a directory structure in README.md - not positive this is what we want yet... 


### Building DB

/project/fas/jetz/data/MPYC/Anthropause/Inputs/MoseyDB

### Process and Clean Data

* Using a control file approach for the study period segmentation/annotation b/c this allows future flexibility to assign individual-/cohort-specific study periods according to some covariate (e.g., taxa, location, etc.)  
* Script complete - saves two new tables to the db in "/processed_data":  
  * *event_mod* has occurences outside stufy periods removed and the basic study periods annotated  
  * *event_clean* further modifies that to remove outliers.  Points >95% quantile for _both_ step length and turn angle are removed.
  
  With the test dataset (~2.5 million occurrences) it's pretty fast (< 2 min) - may need to run on high memory HPC for full dataset


### Analysis
* Potential mods to analysis framework:  
  * Take the approach of measuring the delta in the diffs between periods.  So w/i a year measure diff of X btw pre- and post-ld; then measure change in that diff between 2019 (control) and 2020 (treatment). Can do this for:  
    * SL and TA distributions  
    * UDs - see below, switch to dBBMM???  
    * Niches - breadth and dissimilarity (the latter already being the diff metric)  
  * For HR/UD - I actually think the time-invariant, stationary HR from an AKDE is not really what we want...  also it carries additional assumptions.  But if we fit a dBBMM we get the UD, and can use Earth Movers Distance to compare between UDs.  
    * Not positive but I think we can also fit a single dBBMM per year and chop up by ld period after the fact.  
  * For niches - use spider plots to show relative contribution ot breadth and dissimilarity, and colored traces to show e.g., ld and per-ld

###Annotations
Pick up by running import and annotation steps in wf-mosey_env_BG.sh (after generate bg is complete on HPC)