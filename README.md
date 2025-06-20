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

Data are stored as a [mosey_db](https://github.com/benscarlson/mosey_db), which is a SQLite relational database built to store data from [MoveBank](www.movebank.org). Raw wildlife movement data associated with this project includes data for species that pose conservation concerns and therefore cannot be released publicly. Aggregated data products at the individual-year-week scale, model files and plots, and tabular model results are publicly available [here on Open Science Framework](https://osf.io/3ua2c/files/osfstorage).

## Citation

This project is led by:

Walter Jetz (PI)<sup>1</sup>  
Ruth Oliver<sup>1</sup>  
Diego Ellis-Soto<sup>1</sup>  
Scott Yanco<sup>1</sup>  
Brett Jesmer<sup>1, 2</sup>  

<sup>1</sup>Center for Biodiversity and Global Change  
Department of Ecology and Evolutionary Biology  
Yale University  
New Haven, CT, USA  

<sup>2</sup>Department of Fish and Wildlife Conservation  
Virginia Tech   
Blacksburg, VA, USA  

## Additional Information & Resources

Published comment article describing the research opportunity: [https://www.nature.com/articles/s41559-020-1237-z](https://www.nature.com/articles/s41559-020-1237-z (open access))
