####################################################################################
## Purpose: Process spanish datasets for analysis
##
## Author: Nick Brazeau
##
## Date: 26 May, 2020
####################################################################################
#remotes::install_github("mrc-ide/drjacoby", ref = "develop")
library(drjacoby)
#remotes::install_github("mrc-ide/COVIDCurve")
library(COVIDCurve)

#............................................................
# AGE-Band Analysis
#...........................................................
ageband <- readRDS("data/derived/ESP/ESP_agebands.RDS")

#......................
# get data in format for COVIDCurve
#......................
dat <- list(obs_deaths = ageband$deaths,
            obs_serologyrate = ageband$seroprev$seroprev)

#..................
# make model
#..................
deathdf_params <- data.frame(name = paste0("ma", 1:10),
                             min =  rep(0, 10),
                             init = rep(0.5, 10),
                             max =  rep(1, 10),
                             dsc1 = rep(0, 10),
                             dsc2 = rep(1, 10))

infxndf_params <- tibble::tibble(name = paste0("y", 1:5),
                                 min  = c(0,   0,   0,  0,   0),
                                 init = c(0.5, 0.5, 9,  0.5, 0.5),
                                 max =  c(1,   1,   12, 1,   1),
                                 dsc1 = c(0,   0,    0,  0,   0),
                                 dsc2 = c(1,   1,    12, 1,   1))
knotdf_params <- tibble::tibble(name = paste0("x", 1:4),
                                min  = c(0,    0.30, 0.50, 100),
                                init = c(0.05, 0.50, 0.60, 120),
                                max =  c(0.5,  0.75, 0.99, 137),
                                dsc1 = c(0,    0.30, 0.50, 100),
                                dsc2 = c(0.5,  0.75, 0.99, 137))

# infxndf_params <- data.frame(name = paste0("y", 1:5),
#                              min = c(0,     0,   0,   0,   0),
#                              init = c(0.01, 0.5, 2, 0.5, 0.5),
#                              max = c(0.05,  1,   10, 1,   1),
#                              dsc1 = c(0,    0,   0,   0,   0),
#                              dsc2 = c(0.05, 1,   10, 1,   1))
#
# knotdf_params <- tibble::tibble(name = paste0("x", 1:4),
#                                 min  = c(0,    0.50, 0.80, 100),
#                                 init = c(0.05, 0.60, 0.85, 120),
#                                 max =  c(0.25, 0.90, 0.99, 137),
#                                 dsc1 = c(0,    0.50, 0.80, 100),
#                                 dsc2 = c(0.25, 0.90, 0.99, 137))

serodf_params <- tibble::tibble(name = c("sens", "spec", "sero_rate", "sero_day"),
                                min =  c(0.82,    0.98,   10,         117),
                                init = c(0.85,    0.99,   10,         125),
                                max =  c(0.88,    1,  10,         131),
                                dsc1 = c(8500,     9900,    10,         117),
                                dsc2 = c(1500,     100,     10,         131))

df_params <- rbind.data.frame(deathdf_params, infxndf_params, knotdf_params, serodf_params)

mod1 <- COVIDCurve::make_IFRmodel_agg$new()
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_level("Time-Series")
mod1$set_data(dat)
mod1$set_IFRparams(paste0("ma", 1:10))
mod1$set_maxMa("ma10")
mod1$set_Infxnparams(paste0("y", 1:5))
mod1$set_relInfxn("y3")
mod1$set_Knotparams(paste0("x", 1:4))
mod1$set_relKnot("x4")
mod1$set_Seroparams(c("sens", "spec", "sero_rate", "sero_day"))
mod1$set_popN(ageband$popN)
mod1$set_paramdf(df_params)
mod1$set_pa(ageband$pa$pop_prop)
mod1$set_rcensor_day(.Machine$integer.max) # serostudy may have gone up to 131, data only goes up to 137
#..................
# run model
#..................
start <- Sys.time()
r_mcmc_out.ageband <- COVIDCurve::run_IFRmodel_agg(IFRmodel = mod1,
                                                   reparamIFR = TRUE,
                                                   reparamInfxn = TRUE,
                                                   reparamKnots = TRUE,
                                                   chains = 10,
                                                   burnin = 1e4, samples = 1e3)
Sys.time() - start
# out
saveRDS(r_mcmc_out.ageband, "data/derived/ESP/ESP_mcmc_agebands_fit.rds")

