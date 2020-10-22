#........................................................................................
## Purpose: Use Drake to fit IFR Models
##
## Notes:
#........................................................................................
setwd("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis")
library(drake)
library(parallel)
library(COVIDCurve)
library(tidyverse)
source("R/covidcurve_helper_functions.R")
#...................................................................................
# Make Paramset and write to disk for input into MCMC
#.................................................................................
#......................
# onset to deaths
#......................
tod_paramsdf <- tibble::tibble(name = c("mod", "sod", "sero_con_rate"),
                               min  = c(18,     0,     16),
                               init = c(19,     0.85,  18),
                               max =  c(20,     1,     21),
                               dsc1 = c(19.8,   2550,  18.3),
                               dsc2 = c(0.1,    450,   0.1))

#............................................................
#---- New York City #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.85,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(204.5,   287.5),
                                dsc2 =  c(30.5,    0.5))

#......................
# get fits from stan model
#......................
nyssens <- readr::read_csv("results/prior_inputs/NYS1_sens_reg_age.csv")
sens <- fitdistrplus::fitdist(unlist(nyssens), distr = "beta", method = "mme")
nysspec <- readr::read_csv("results/prior_inputs/NYS1_spec_reg_age.csv")
spec <- fitdistrplus::fitdist(unlist(nysspec), distr = "beta", method = "mme")
sens_spec_tbl$dsc1[sens_spec_tbl$name == "sens"] <- sens$estimate[["shape1"]]
sens_spec_tbl$dsc2[sens_spec_tbl$name == "sens"] <- sens$estimate[["shape2"]]
sens_spec_tbl$dsc1[sens_spec_tbl$name == "spec"] <- spec$estimate[["shape1"]]
sens_spec_tbl$dsc2[sens_spec_tbl$name == "spec"] <- spec$estimate[["shape2"]]


#......................
# agebands
#......................
rawage <- readRDS("data/derived/USA/NYC1_agebands.RDS")
NYC_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 5, maxMa = "ma5",
                                            groupvar = "ageband",  dat = rawage,
                                            num_xs = 4, max_xveclist = list("name" = "x4", min = 214, init = 223, max = 230, dsc1 = 214, dsc2 = 230),
                                            num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 15.94, dsc1 = 0, dsc2 = 15.94),
                                            sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)




#............................................................
# MCMC Run
#...........................................................
bvec <- seq(5, 2.5, length.out = 50)

fit <- COVIDCurve::run_IFRmodel_age(IFRmodel = NYC_age_mod,
                                    reparamIFR = TRUE,
                                    reparamInfxn = TRUE,
                                    reparamKnots = TRUE,
                                    binomial_likelihood = TRUE,
                                    chains = 10,
                                    burnin = 1e4,
                                    samples = 1e4,
                                    rungs = 50,
                                    GTI_pow = bvec,
                                    cluster = cl,
                                    thinning = 10)

# save out
dir.create("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/NYC_test/", recursive = TRUE)
outpath <- "/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/results/NYC_test/NYC1_noserorev_fit.RDS"
saveRDS(fit, file = outpath)
