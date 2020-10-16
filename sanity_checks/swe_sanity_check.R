#........................................................................................
## Purpose: sweden definitive death sanity check
##
## Notes:
#........................................................................................
setwd("/proj/ideel/meshnick/users/NickB/Projects/neil_def_deaths/reestimate_covidIFR_analysis")
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
#---- SWE1 #----
#...........................................................
sens_spec_tbl <- tibble::tibble(name =  c("sens", "spec"),
                                min =   c(0.50,    0.50),
                                init =  c(0.99,    0.99),
                                max =   c(1.00,    1.00),
                                dsc1 =  c(156.5,   267.5),
                                dsc2 =  c(1.5,     3.5))
#......................
# agebands
#......................
rawage <- readRDS("data/derived/SWE1/SWE1_agebands.RDS")
SWE_age_mod <- make_noSeroRev_IFR_model_fit(num_mas = 10, maxMa = "ma10",
                                            groupvar = "ageband",  dat = rawage,
                                            num_xs = 4, max_xveclist = list("name" = "x4", min = 214, init = 223, max = 230, dsc1 = 214, dsc2 = 230),
                                            num_ys = 5, max_yveclist = list("name" = "y3", min = 0, init = 9, max = 16.15, dsc1 = 0, dsc2 = 16.15),
                                            sens_spec_tbl = sens_spec_tbl, tod_paramsdf = tod_paramsdf)


#............................................................
# run mcmc
# going to do just casual 1e4/1e4 iterations
#...........................................................
n_chains <- 10
n_cores <- parallel::detectCores()

if (n_cores < n_chains) {
  mkcores <- n_cores - 1
} else {
  mkcores <- n_chains
}

cl <- parallel::makeCluster(mkcores)

# logit cases
fit <- COVIDCurve::run_IFRmodel_age(IFRmodel = SWE_age_mod,
                                    reparamIFR = TRUE,
                                    reparamInfxn = TRUE,
                                    reparamKnots = TRUE,
                                    binomial_likelihood = FALSE,
                                    chains = n_chains,
                                    burnin = 1e4,
                                    samples = 1e4,
                                    rungs = 50,
                                    GTI_pow = list(seq(5, 2.5, length.out = 50)),
                                    cluster = cl,
                                    thinning = 10)


# out
dir.create("/proj/ideel/meshnick/users/NickB/Projects/neil_def_deaths/reestimate_covidIFR_analysis/results/Modfits_noserorev/", recursive = TRUE)
outpath = paste0("/proj/ideel/meshnick/users/NickB/Projects/neil_def_deaths/reestimate_covidIFR_analysis/results/Modfits_noserorev/",
                 mod$name, "_rung", mod$rungs, "_burn", mod$burnin, "_smpl", mod$samples, "_NoSeroRev.RDS")
saveRDS(fit, file = outpath)
