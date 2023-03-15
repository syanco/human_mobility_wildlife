## TODO

-   Export FINAL conda env .yml files
    -   covid\
    -   spatial\
-   Clean comments on script 01\
-   Check that dBBMM log file count matches the out/ directory

## Notes

\*Fitting dBBMMs for each ind x yr combination.

\*I'm getting an error in thr prep and clean script on the HPC - that I DO NOT GET when run local:

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

For now I'm just scp'ing the processed db to the HPC to run the dbbmms

-   For the HR size - consider just a lmm of size \~ trt X yr + (1\|ind) as basic model

### General

-   Proposed a directory structure in README.md - not positive this is what we want yet...

### Building DB

/project/fas/jetz/data/MPYC/Anthropause/Inputs/MoseyDB

### Process and Clean Data

-   Using a control file approach for the study period segmentation/annotation b/c this allows future flexibility to assign individual-/cohort-specific study periods according to some covariate (e.g., taxa, location, etc.)\

-   Script complete - saves two new tables to the db in "/processed_data":

    -   *event_mod* has occurences outside stufy periods removed and the basic study periods annotated\
    -   *event_clean* further modifies that to remove outliers. Points \>95% quantile for *both* step length and turn angle are removed.

    With the test dataset (\~2.5 million occurrences) it's pretty fast (\< 2 min) - may need to run on high memory HPC for full dataset

### Analysis

-   Potential mods to analysis framework:
    -   Take the approach of measuring the delta in the diffs between periods. So w/i a year measure diff of X btw pre- and post-ld; then measure change in that diff between 2019 (control) and 2020 (treatment). Can do this for:
        -   SL and TA distributions\
        -   UDs - see below, switch to dBBMM???\
        -   Niches - breadth and dissimilarity (the latter already being the diff metric)\
    -   For HR/UD - I actually think the time-invariant, stationary HR from an AKDE is not really what we want... also it carries additional assumptions. But if we fit a dBBMM we get the UD, and can use Earth Movers Distance to compare between UDs.
        -   Not positive but I think we can also fit a single dBBMM per year and chop up by ld period after the fact.\
    -   For niches - use spider plots to show relative contribution ot breadth and dissimilarity, and colored traces to show e.g., ld and per-ld

###Annotations Pick up by running import and annotation steps in wf-mosey_env_BG.sh (after generate bg is complete on HPC)

## Workflow Log

Log describes Feb 2023 re-run of the workflow after update to the mosey.db.

### 2023 February 10

Need to export new dataset to Ben and Fede.

Started workflow:\
- moved raw db to processed dir.\
- Ran cleaning script

### 2023 March 1

Got the import tweaked and completed for annotating the mosey swap.

### 2023 March 2

Got the next step done (ran annos for tmax and ndvi) (had to make soem tweaks to make mosey work - basically was able to template from my dynamic niches repo - primar issue was that the current build of mosey tries to loop over entities in both the outer .sh script and the inner r script. I killed the loop in the r script b/c it makes nice print statements basically)

### 2023 March 3

Fixing up import script just like the last two steps - gotta copy changes to dynamic niches. I forgot thought to fix the thing where the anno is saving each csv to int's own filter. better to have fixed at that step but it takes a day to annotate so I don't want to re-run.

OK - got it running. Next step is to grab the SG annotations and then do the merge.

About to walk away for the day and not sure the loop is going to re-start, so the last ind I see printed is.... 2561430505 2561486584

### 2023 March 7

Was able to finally complete the mosey re-annos on the swap db yesterday. Starting the SG annotation steps today - so far seems to work fin, just need to update the filepaths to db (call the swap db and also move off of Loomis filesystem)

Workflow run through the annotate-events-ghm step.

### 2023 March 8

Working on merge script.

Since I'm going to overwrite tables in the main mosey - I made a copy in processed data (mosey_mod_copy.db)

Merge ran successfully on the HPC - back to re-running analyses which **should** run out of the box....

### 2023 March 10

Re-running analysis, completed today: - fit dBBMMs\
- calc space use

Had to restart the fit dbmms b/c the log was saving out the ind_id as scientific notation and it was breaking all the ctf lookups downstream (e.g. in calc space use)

Calc space use is still coming up empty - not sure what's going on - need to debug next week...

### 2023 March 13

Start trying to debug why calc space use is coming up empty...

Fit dbbmm was not upated correctly - added `as.integer()` wrapper to all the outlogs (was only doing it for those with too few events.)\
Relaunched fit_dbbmm.r at 10am

That run didn't work - but did some more debugging - i think the db merge was creating duplicates in some of the tables? Either way I was getting duplicated results, so I added some checks for duoplicates in the fitdbbmm script and modified how I write out to the log. Launcehd a new run at 1500.

### 2023 March 14

Started a new run of calc size to see if yesterdays changes helped.... not sure they did...

I really think this comes back to the db. When I log in and do an OOD session I can't even connect to it - I get malformed event_cbg already exists... (not sure why that doesn't cause an error when run via SLURM though...)

I really think I need to restart this whole merge thing again. This tiume I would mergee upfront and reannotate all in one go. If I do it that way i need to add study ID to diego's swap db...

### 2023 March 14

Starting the merge from scratch, there is a bit of confusion about which db is the right one to start with....

Here's what I'm doing:\
\* Pulled mosey.db from the from_loomis/ dir on HPC to raw_data/\
\* Copied it to processed_data/mosey_mod.db, next to mosey_swap.db from Diego.\
\* Run poc/merge_dbs.r (script that merges the swap db into the study, individual, and event tables of the main mosey.)

