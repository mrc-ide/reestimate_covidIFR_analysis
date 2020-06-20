# Reestimating IFRs with Serology Data Research Compendium 

The structure of this research compendium loosely follows the framework of an R-package for reproducing our analysis in line with the framework described by the [ROpenSci Team](https://github.com/ropensci/rrrpkg).  

An overview of the compendium is below: 
```
├── R                                 # "Helper" R functions specific to this research project 
│   ├── [a-zA-Z0-9_.-].R
├── README.md
├── analyses                            # Scripts for data wrangling, fitting, figures, and tables
│   ├── descriptive_data_plots.R           # Generate descriptive plots
│   └── run_process_country_data.R         # Process data for model fitting
├── data                                 # DO NOT EDIT ANY FILES IN THIS DIRECTORY BY HAND
│   ├── raw/
│   ├── derived
├── figures                              # Figures for Manuscript (and extras)
│   ├── [a-zA-Z_].jpg
├── drake_make.R                         # Drake "Make" Script to generate project 
├── reestimate_covidIFR_analysis.Rproj   # R-proj
├── reports                              # Rmd files for Country-Specific Fits
│   ├── [A-Z]_IFRcountry_results.Rmd        # Individual "sub"-reports
│   ├── main.Rmd                            # Main report that knits together "sub"-reports
│   ├── sim_pwr.Rmd
├── results                              # Results from slurm runs
├── slurm_workers                        # Model fits for running on a slurm cluster
│   ├── run_[A-za-z].R                      # Individual model fits by location
└── tests                                # Temporary developer tests
    └── test_drake.R
```
