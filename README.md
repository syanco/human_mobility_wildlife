# COVID-19 Lockdowns & Animal Movement

Repository containing code for analysis of animal movement in relation to COVID-19 lockdowns in the United States. This is part of the larger "*COVID-19 Bio-logging Initiative*", organized by the International [Bio-logging Initiative](www.bio-logging.net).

Contact information for this repo:  yancos@si.edu

## Code Usage

The code in this repository anticipates the following minimal directory structure to run:

```
/human_moblity_wildlife  
  |  
  +--/src               # Source code directory
  |   |   
  |   +--/funs          # Source code for custom functions called by other scripts
  |   |   
  |`  +--/hpc           # Scripts for submitting jobs to Slurm manager 
  |   | 
  |   +--/workflow      # Scripts that execute elements of the workflow
  |   |
  |   +--/mosey         # Scripts for database annotation
  |   |
  |   +--/startup.R     # Source code for some basic environment configuration
  |    
  +--/ctfs              # Control files  
  | 
  +--/conda_envs        # Stores .yml files with conda environment specifications
  |  
  +--/raw_data          # Raw data stored as initially received, inlcuding database
  |  
  +--/processed_data    # Processed data products, including the working version of the database
  |  
  +--/out               # Analytical outputs, interim products
  |  
  +--/figs              # Scripts to produce figures and outputs

```

Data are stored as a [mosey_db](https://github.com/benscarlson/mosey_db), which is a SQLite relational database built to store data from [MoveBank](www.movebank.org).

[This repository release](https://github.com/julietcohen/mosey_db/releases/tag/v1.0.0) includes the forked `mosey_db` code used to build the database used as input for this repository's workflow.

Raw wildlife movement data associated with this project includes data for species that pose conservation concerns and therefore cannot be released publicly. Aggregated data products at the individual-year-week scale, model files and plots, and tabular model results are publicly available [here on Open Science Framework](https://osf.io/3ua2c/files/osfstorage).

## Additional Information & Resources

Published comment article describing the research opportunity: [https://www.nature.com/articles/s41559-020-1237-z](https://www.nature.com/articles/s41559-020-1237-z (open access))
