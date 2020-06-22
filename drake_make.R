##....................................................................................
## Purpose: Make file for re-estimating IFRs with serology updates
##
## Author: Nick Brazeau
##
## Date: 07 June, 2020
##....................................................................................
library(drake)
source("R/assertions_v5.R")
source("R/process_data.R")
library(drjacoby)
library(COVIDCurve)

#### Setup #####
dir.create("out/tables", recursive = T)
dir.create("out/tables", recursive = T)

#### Make Drake Plan ####
plan <- drake::drake_plan(
  #............................................................
  # assertions that files exist
  #...........................................................

  #............................................................
  # data input for Manuscript
  #...........................................................
  #......................
  # Spain
  #......................
  ESP.agebands.dat = process_data(deaths = "https://www.dropbox.com/s/z42zfiyr9tr86qe/deaths_v2.csv?dl=1",
                                  population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                  sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                  seroprev = "https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1",
                                  cumulative = TRUE,
                                  ECDC = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                  groupingvar = "ageband",
                                  study_ids = "ESP1",
                                  ecdc_countrycode = "ESP",
                                  filtRegions = NULL,
                                  filtGender = NULL,
                                  filtAgeBand = c("0-10", "10-20", "20-30",
                                                  "30-40", "40-50", "50-60",
                                                  "60-70", "70-80", "80-90", "90-999")),

  #............................................................
  # Param Df for Manuscript Fits
  #...........................................................
  #......................
  # read in paramdf and split by fit
  #......................


  #............................................................
  # plots & tables for manuscript
  #...........................................................


  #............................................................
  # Simulation Studies for Supplement
  #...........................................................

  #............................................................
  # render supplement
  #...........................................................
  report = rmarkdown::render(
    knitr_in("Supplement/SupplementaryMaterials.Rmd"),
    output_file = file_out("SupplementaryMaterials.Rmd"),
    quiet = TRUE
  )
)


#### Viz Drake Plan ####
ortrta
drake::vis_drake_graph(ortrta)

#### Run Drake Plan ####
drake::make(ortrta)


