####################################################################################
## Purpose: Process spanish datasets for analysis
##
## Author: Nick Brazeau
##
## Date: 26 May, 2020
####################################################################################
#remotes::install_github("mrc-ide/drjacoby", ref = "develop")
library(drjacoby)
#remotes::install_github("mrc-ide/COVIDCurve", ref = "seroconv_ext")
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
                             min =  rep(0, 10), rep(0, 10),
                             init = rep(0.5, 10), rep(7, 10),
                             max =  rep(1, 10), rep(15, 10),
                             dsc1 = rep(0, 20),
                             dsc2 = rep(1, 10), rep(15, 10))


serodf_params <- tibble::tibble(name = c("sens", "spec", "sero_rate", "sero_date"),
                                min =  c(0.82,    0.97,   0.1,         117),
                                max =  c(0.88,    0.99,   0.1,         131),
                                init = c(0.85,    0.98,   0.1,         125),
                                dsc1 = c(850,     980,    100,         125),
                                dsc2 = c(150,     20,     900,         5))

df_params <- rbind.data.frame(deathdf_params, serodf_params)

mod1 <- make_modinf_agg$new()
mod1$set_level("Time-Series")
mod1$set_data(dat)
mod1$set_IFRparams(c(paste0("r", 1:9), "ma10"))
mod1$set_maxMa("ma10")
mod1$set_Infxnparams(c("y1", "y2", "y3"))
mod1$set_Seroparams(c("sens", "spec", "sero_rate", "sero_date"))
mod1$set_popN(sum(ageband$popN))
mod1$set_paramdf(df_params)
mod1$set_pa(ageband$pa$pop_prop)
mod1$set_MeanOnset(18.8)
mod1$set_CoefVarOnset(0.45)
mod1$set_knots(c(1, 90, 132))

#..................
# run model
#..................
r_mcmc_out.ageband <- COVIDCurve::run_modinf_agg(modinf = mod1, reparamIFR = TRUE, rungs = 10)
# out
saveRDS(r_mcmc_out.ageband, "data/derived/ESP/ESP_mcmc_agebands.rds")
saveRDS(mod1, "data/derived/ESP/ESP_modinf_agebands.rds")

