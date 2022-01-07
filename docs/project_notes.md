## TODO  
* Export conda env .yml files

## Major Activity Log

|Date|Activity|
|:-|:------------|
|2022-01-05|Initialize repo, docs, etc.|
|2022-01-07|Start building top-level workflow and cleaning scripts|


## Notes  

Left off at outlier removal in script 01 - pick up there.

### General

* Proposed a directory structure in README.md - not positive this is what we want yet... 


### Building DB


### Process and Clean Data

* Using a control file approach for the study period segmentation/annotation b/c this allows future flexibility to assign individual-/cohort-specific study periods according to some covariate (e.g., taxa, location, etc.)


### Analysis
* Potential mods to analysis framework:  
  * Take the approach of measuring the delta in the diffs between periods.  So w/i a year measure diff of X btw pre- and post-ld; then measure change in that diff between 2019 (control) and 2020 (treatment). Can do this for:  
    * SL and TA distributions  
    * UDs - see below, switch to dBBMM???  
    * Niches - breadth and dissimilarity (the latter already being the diff metric)  
  * For HR/UD - I actually think the time-invariant, stationary HR from an AKDE is not really what we want...  also it carries additional assumptions.  But if we fit a dBBMM we get the UD, and can use Earth Movers Distance to compare between UDs.  
    * Not positive but I think we can also fit a single dBBMM per year and chop up by ld period after the fact.  
  * For niches - use spider plots to show relative contribution ot breadth and dissimilarity, and colored traces to show e.g., ld and per-ld
    