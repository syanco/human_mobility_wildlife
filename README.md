# COVID-19 Lockdowns & Animal Movement

Repository containing code for analysis of animal movement in relation to COVID-19 lockdowns in the United States. This is part of the larger "*COVID-19 Bio-logging Initiative*", organized by the International Bio-loggin Initiative ([www.bio-logging.net](www.bio-logging.net)).

Contact information for this repo:  scott.yanco@yale.edu

## Code Usage

The code in this repository anticipates the following minimal directory structure to run:

```
/project_root  
  |  
  +--/analysis          # *This repository*  
  |   |  
  |   +--/docs          # Supporting documentation (other than README.md)  
  |   |  
  |   +--/src           # Source code directory, includes top level workflow script
  |   |   |
  |   |   +--/funs      # Source code for custom functions called by other scripts
  |   |   | 
  |`  |   +--/hpc       # Scripts for submitting jobs to Slurm manager
  |   |   |
  |   |   +--/workflow  # Scripts that execute elements of the workflow
  |   |   |
  |   |   +--/poc       # Scripts that are in active development or scratch scripts excl;uded from workflow
  |   |   |
  |   |   ...           # Other directories as needed for particular non-workflow tasts (e.g., method validation analyses, interim vizualization scripts, etc.) 
  |   |  
  |   +--/ctfs          # Control files  
  |   |
  |   +--/conda_envs    # Stores .yml files wiwth conda environment recipes
  |  
  +--/raw_data          # Raw data stored as initially received, inlcuding mosey_db
  |  
  +--/processed_data    # Any processed data products, especially the wokring version of the mosey db
  |  
  +--/out               # Analytical outputs, interim products
  |  
  +--/figs              # Figures produced by this workflow
  |  
  ...                   # Other scripts as needed
```

Data associated with this project are far too large to be stored on GitHub and are currently private anyhow.  Raw data are stored as a [mosey_db](https://github.com/benscarlson/mosey_db), which is a SQLite relational database built to store data from [MoveBank](www.movebank.org). If data become publically available the archive will be described here (or the code to build the mosey_db from MoveBank will be released).

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

Recently published comment article describing the opportunity: [https://www.nature.com/articles/s41559-020-1237-z](https://www.nature.com/articles/s41559-020-1237-z (open access))
