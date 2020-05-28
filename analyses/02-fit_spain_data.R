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
deathdf_params <- data.frame(name = c(paste0("r", 1:9), "ma10", paste0("y", 1:10)),
                             min =  c(rep(0, 10), rep(0, 10)),
                             init = c(rep(0.5, 10), rep(7, 10)),
                             max =  c(rep(1, 10), rep(17, 10)),
                             dsc1 = c(rep(0, 20)),
                             dsc2 = c(rep(1, 10), rep(17, 10))
                             )


serodf_params <- tibble::tibble(name = c("sens", "spec", "sero_rate", "sero_date"),
                                min =  c(0.82,    0.98,   10,         117),
                                max =  c(0.88,    0.999,  10,         131),
                                init = c(0.85,    0.99,   10,         125),
                                dsc1 = c(850,     990,    10,         117),
                                dsc2 = c(150,     10,     10,         131))

df_params <- rbind.data.frame(deathdf_params, serodf_params)

mod1 <- make_modinf_agg$new()
mod1$set_level("Time-Series")
mod1$set_data(dat)
mod1$set_IFRparams(c(paste0("r", 1:9), "ma10"))
mod1$set_maxMa("ma10")
mod1$set_Infxnparams(paste0("y", 1:10))
mod1$set_Seroparams(c("sens", "spec", "sero_rate", "sero_date"))
mod1$set_popN(ageband$popN)
mod1$set_paramdf(df_params)
mod1$set_pa(ageband$pa$pop_prop)
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_knots(c(1, 16, 30, 45, 59, 74, 88, 103, 117, 132))

#..................
# run model
#..................
start <- Sys.time()
r_mcmc_out.ageband <- COVIDCurve::run_modinf_agg(modinf = mod1, reparamIFR = TRUE,
                                                 chains = 5, rungs = 25,
                                                 GTI_pow = 3.0,
                                                 burnin = 1e3, samples = 1e3)
Sys.time() - start
# out
saveRDS(r_mcmc_out.ageband, "data/derived/ESP/ESP_mcmc_agebands_25rungs_5chains_3GTI.rds")
saveRDS(mod1, "data/derived/ESP/ESP_modinf_agebands.rds")

