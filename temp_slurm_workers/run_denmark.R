####################################################################################
## Purpose: run Denmark on LL
##
## Notes:
####################################################################################
setwd("/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis")
library(rslurm)
library(COVIDCurve)
library(tidyverse)
source("R/covidcurve_helper_functions.R")

#............................................................
# read in AGE DATA
#...........................................................
age <- readRDS("data/derived/DNK/DNK_agebands.RDS")

#......................
# make agemod for age
#......................
ifr_paramsdf <- make_ma_reparamdf(num_mas = 5)
knot_paramsdf <- make_splinex_reparamdf(max_xvec = list("name" = "x4", min = 125, init = 131, max = 138, dsc1 = 125, dsc2 = 138),
                                        num_xs = 4)
infxn_paramsdf <- make_spliney_reparamdf(max_yvec = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                                         num_ys = 5)
sero_paramsdf <- tibble::tibble(name =  c("sens", "spec", "sero_rate", "sero_day"),
                                min =   c(0.80,     0.97,   10,         97),
                                init =  c(0.82,     0.99,   10,         100),
                                max =   c(0.84,     1.00,   10,         105),
                                dsc1 =  c(8200,     9900,    5,         100),
                                dsc2 =  c(1800,     100,     15,        0.1))

df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sero_paramsdf)

# format data
inputdata <- list(obs_deaths = age$deaths,
                  obs_serologyrate = age$seroprev$seroprev)
# make pa
pa <- age$prop_pop$pop_prop

#......................
# make model
#......................
agemod <- COVIDCurve::make_IFRmodel_agg$new()
agemod$set_MeanOnset(18.8)
agemod$set_CoefVarOnset(0.45)
agemod$set_level("Time-Series")
agemod$set_data(inputdata)
agemod$set_IFRparams(paste0("ma", 1:5))
agemod$set_maxMa("ma5")
agemod$set_Knotparams(paste0("x", 1:4))
agemod$set_relKnot("x4")
agemod$set_Infxnparams(paste0("y", 1:5))
agemod$set_relInfxn("y3")
agemod$set_Seroparams(c("sens", "spec", "sero_rate", "sero_day"))
agemod$set_popN(age$popN)
agemod$set_paramdf(df_params)
agemod$set_pa(pa)
agemod$set_rcensor_day(.Machine$integer.max)

#............................................................
# read in REGION DATA
#...........................................................
region <- readRDS("data/derived/DNK/DNK_regions.RDS")

#......................
# make regionmod for region
#......................
ifr_paramsdf <- make_ma_reparamdf(num_mas = 3)
# save out key for regions
rgnkey <- tibble::tibble(rgn = region$deaths$region[region$deaths$ObsDay == 1], name = ifr_paramsdf$name)
saveRDS(rgnkey, "data/derived/DNK/DNK_rgn_key.RDS")

df_params <- rbind.data.frame(ifr_paramsdf, infxn_paramsdf, knot_paramsdf, sero_paramsdf)

# format data
inputdata <- list(obs_deaths = region$deaths,
                  obs_serologyrate = region$seroprev$seroprev)
# make pa
pa <- region$seroprev_group$seroprev
pa <- pa/sum(pa)
#......................
# make model
#......................
regionmod <- COVIDCurve::make_IFRmodel_agg$new()
regionmod$set_MeanOnset(18.8)
regionmod$set_CoefVarOnset(0.45)
regionmod$set_level("Time-Series")
regionmod$set_data(inputdata)
regionmod$set_IFRparams(paste0("ma", 1:3))
regionmod$set_maxMa("ma1")
regionmod$set_Knotparams(paste0("x", 1:4))
regionmod$set_relKnot("x4")
regionmod$set_Infxnparams(paste0("y", 1:5))
regionmod$set_relInfxn("y3")
regionmod$set_Seroparams(c("sens", "spec", "sero_rate", "sero_day"))
regionmod$set_popN(age$popN)
regionmod$set_paramdf(df_params)
regionmod$set_pa(pa)
regionmod$set_rcensor_day(.Machine$integer.max)


#......................
# simple wrapper for run on slurm
#......................
run_wrapper <- function(modelobj, GTI_pow) {
  COVIDCurve::run_IFRmodel_agg(IFRmodel = modelobj,
                               reparamIFR = TRUE,
                               reparamInfxn = TRUE,
                               reparamKnots = TRUE,
                               chains = 10,
                               burnin = 1e4,
                               samples = 1e4,
                               rungs = 25,
                               GTI_pow = GTI_pow)
}
#......................
# make "map" for slurm
#......................
map <- tibble::tibble(modelobj = list(agemod, regionmod), GTI_pow = 3)

#......................
# send out on slurm
#......................
ntry <- 3
sjob <- rslurm::slurm_apply(f = run_wrapper,
                            params = map,
                            jobname = 'DNKfits',
                            nodes = ntry,
                            cpus_per_node = 1,
                            submit = T,
                            slurm_options = list(mem = "24g",
                                                 'cpus-per-task' = 1,
                                                 error =  "%A_%a.err",
                                                 output = "%A_%a.out",
                                                 time = "1-00:00:00"))

















# sanity
