## TODO

-   Export FINAL conda env .yml files
    -   covid\
    -   spatial\
-   Clean comments on script 01\
-   Check that dBBMM log file count matches the out/ directory

## Current Working Notes



### General

-   Proposed a directory structure in README.md - not positive this is what we want yet...


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

"Final" disposition of dbs...

**raw_data/**
mosey.db - the pre-swap mosey (from April 2022)
mosey_swap.db - the db with replacement studies that Diego made.
mosey_merge.db - the merged, complete raw db.

**processed_data/**
mosey_mod.db - the working version of the complete db, by tomorrow, everything should run form this, typically using the 'event_clean' table
mosey_swap.db - copy of above, used as working version


### 2023 March 17

Successfully re-ran clean script


### 2023 March 20

Starting annotation again - Ts & Ps for me.

Very minor edits needed to get the ingest step working (switched to using study_id as entity rather than individual)

### 2023 April 26

Finished several key re-runs over the last few days:

1). *Main model selection scripts*

Using the following approach:

We a priori expected an interaction effect - so we fit an interaction model for every species.  In cases where the interaction effect was absent (credible intervals overlapping zero) we re-fit the model considering additive effects of SG and GHM only and present the results of only that model.

Successful re-fit models with the new build of the db, and collated results - passed on for viz.

2). *Intra individual models*

Added sapce use as control var in the niche models. cleanded up HPC scripts and relaunched on HPC (currenlty running)

