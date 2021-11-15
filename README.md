# Reestimating IFRs with Serology Data Research Compendium 
[![DOI](https://zenodo.org/badge/267118325.svg)](https://zenodo.org/badge/latestdoi/267118325)

The structure of this research compendium loosely follows the framework of an R-package for reproducing our analysis (for more information on research compendiums, see [ROpenSci Team](https://github.com/ropensci/rrrpkg)).  

An overview of the compendium is below: 
```
├── R                                 # "Helper" R functions specific to this research project 
├── README.md
├── analyses/                         # Scripts for data wrangling, fitting, figures, and tables
│   ├── ModFits/
│   │   ├── CareHomes_drakeworker.R
│   │   ├── ModFits_drakeworker_v2.R
│   │   ├── seroRev_ModFits_drakeworker_v2.R
│   ├── Rgn_Mod_Stan/
│   │   ├── run_regional_model.R
│   ├── SimWork/
│   │   ├── SeroRev_SimCurves_drakeworker_v2.R
│   │   ├── SimCurves_drakeworker_v2.R
│   │   ├── fit_conceptual_figure_delayeffects.R
│   │   ├── fit_seroday_concept.R
|
│   ├── 00-convalescent_seroAbs_v2.R
│   ├── 01-run_process_country_data_v2.R
│   ├── 02-descriptive_data_plots_v3.R
│   ├── 03-plot_methods_conceptual.R
│   ├── 04-analyze_overall_IFRs.R
│   ├── 05-analyze_age_IFRs_v2.R
│   ├── 06-analyze_carehome_IFRs.R
│   ├── 07-mod_posteriors.R
│   ├── 08-model_ppcs.R
│   ├── 09-collate_sims.R
│   ├── 10-plot_delay_concept.R
│   ├── 11-plot_seroday_concept.R
|
├── data/                                 # DO NOT EDIT ANY FILES IN THIS DIRECTORY BY HAND
│   ├── raw/
│   ├── derived/
│   ├── simdat/
│   ├── plot_aesthetics/
|
├── drake_clst/                          # Instructions for a slurm scheduler for resource allocation
|
├── reports/                             # Reports for simulations and study-specific runs                    
├── reestimate_covidIFR_analysis.Rproj   # R-proj
├── run_nightly_workers.sh               # Workers for generating MCMCMC fits
|                              
```

### Running the Code
Users will first need to perform model fitting for each of the included studies and/or the simulation runs with code provided in the `analyses/ModFits/` and `analyse/SimWork/` directories, respectively. These age-based models depend on the [`COVIDCurve` R-Package](https://github.com/mrc-ide/COVIDCurve). Individual "worker" scripts for models without and with seroconversion are available (you can send these out on a slurm cluster with the `run_nightly_workers.sh` script). From experience on a Linux-based cluster, fits take approximately 8-32 hours to complete. 
  
  
After you have produced the model fits, you can run the rest of the analyses in sequential order 01-11 (script `00-convalescent_seroAbs_v2.R` relies on data that is available upon request. The needed parameter results are provided in the `results/prior_inputs` directory). 


