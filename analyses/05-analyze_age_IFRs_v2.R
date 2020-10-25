##................................................................................................
## Purpose: Analyze age specific IFRs
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")
source("R/covidcurve_helper_functions.R")

# colors now based on location
locatkey <- readr::read_csv("data/plot_aesthetics/color_studyid_map.csv")
mycolors <- locatkey$cols
names(mycolors) <- locatkey$location
# order
order <- readr::read_csv("data/plot_aesthetics/study_id_order.csv")



#............................................................
#---- Read in Fitted and observed data #----
#...........................................................
regrets <- list.files("results/Modfits_noserorev/", full.names = T)
regretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(regrets), "_age", simplify = T)[,1]),
                            sero = "reg",
                            path = regrets)

serorevrets <- list.files("results/ModFits_SeroRev/", full.names = T)
serorevretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(serorevrets), "_age", simplify = T)[,1]),
                                sero = "serorev",
                                path = serorevrets)

# read in observed data
dscdat <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS")
dsc_agedat <- dscdat %>%
  dplyr::filter(breakdown == "ageband") %>%
  dplyr::filter(!grepl("_nch", study_id))

#............................................................
#---- Long Table Age Specific Results #----
#...........................................................
#......................
# get modelled stratified IFRs
#......................
retmapIFR <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(strataIFRret = purrr::map(path, get_strata_IFRs))

modIFR_age <- retmapIFR %>%
  tidyr::unnest(cols = "strataIFRret") %>%
  dplyr::select(c("study_id", "strata", "median", "LCI", "UCI", "sero")) %>%
  dplyr::rename(ageband = strata) %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2))

modIFR_age_column <- modIFR_age %>%
  dplyr::mutate(IFR_cols = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "ageband", "IFR_cols", "sero")) %>%
  tidyr::pivot_wider(., id_cols = c("study_id", "ageband"),
                     names_from = "sero", values_from = "IFR_cols") %>%
  dplyr::rename(regIFR = reg,
                serorevIFR = serorev)

#......................
# get modelled seroprevs
#......................
set.seed(48)
retmapSero <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(strataIFRret = purrr::map(path, get_strata_seroprevs))

#......................
# get simple seropred
#......................
simp_seroprevdat <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "seroprev_adjdat")) %>%
  tidyr::unnest(cols = "seroprev_adjdat") %>%
  dplyr::group_by(study_id, location) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # subset to latest serostudy
  dplyr::filter(!is.na(n_positive) & !is.na(n_tested)) %>%
  dplyr::mutate(sero_midday = seromidpt + lubridate::ymd("2020-01-01") - 1,
                seroprev = n_positive/n_tested) %>%
  dplyr::select(c("study_id", "location", "ageband", "sero_midday", "seroprev"))
# manual liftover for studies that only reported CIs
ci_simp_seroprevdat <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "seroprev_adjdat")) %>%
  tidyr::unnest(cols = "seroprev_adjdat") %>%
  dplyr::group_by(study_id, location) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # subset to latest serostudy
  dplyr::filter(is.na(n_positive) & is.na(n_tested)) %>%
  dplyr::mutate(sero_midday = seromidpt + lubridate::ymd("2020-01-01") - 1) %>%
  dplyr::select(c("study_id", "location", "ageband", "sero_midday", "seroprev"))

# bring together
simp_seroprevdat <- dplyr::bind_rows(simp_seroprevdat, ci_simp_seroprevdat) %>%
  dplyr::rename(obs_input_seropev = seroprev) %>%
  dplyr::mutate(obs_input_seropev = round(obs_input_seropev * 100, 2))


#......................
# now to get inferred seroprevalence
#......................
# am going to need data dictionaries
datdict <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(dictkey = purrr::map(path, get_data_dict)) %>%
  tidyr::unnest(cols = "dictkey")

# am also going to need the seromiddays for comparison
middays <- simp_seroprevdat %>%
  dplyr::select(c("study_id", "location", "ageband", "sero_midday")) %>%
  dplyr::mutate(seromidpt = as.numeric(sero_midday - lubridate::ymd("2020-01-01")) + 1) %>%
  dplyr::select(-c("sero_midday"))

# now out
seroprev_age_column <- retmapSero %>%
  dplyr::select(-c("path")) %>%
  tidyr::unnest(cols = "strataIFRret") %>%
  dplyr::left_join(., middays, by = "study_id") %>%
  dplyr::filter(ObsDay == seromidpt) %>%
  dplyr::select(c("study_id", "sero", "sim", dplyr::starts_with("RG_pd_seroprev_"))) %>%
  tidyr::pivot_longer(., cols = -c("study_id", "sero", "sim"),
                      names_to = "Strata", values_to = "RG_pd_seroprev") %>%
  dplyr::mutate(Strata = gsub("RG_pd_seroprev_", "", Strata)) %>%
  dplyr::left_join(., datdict, by = c("study_id", "Strata", "sero")) %>%
  dplyr::group_by_at(c("study_id", "sero", "ageband")) %>%
  dplyr::summarise(
    LCI = quantile(RG_pd_seroprev, 0.025, na.rm = T),
    median = median(RG_pd_seroprev, na.rm = T),
    UCI = quantile(RG_pd_seroprev, 0.975, na.rm = T))  %>%
  dplyr::mutate(median = round(median * 100, 2),
                LCI = round(LCI * 100, 2),
                UCI = round(UCI * 100, 2)) %>%
  dplyr::mutate(IFR_cols = paste0(median, " (", LCI, ", ", UCI, ")")) %>%
  dplyr::select(c("study_id", "ageband", "IFR_cols", "sero")) %>%
  tidyr::pivot_wider(., id_cols = c("study_id", "ageband"),
                     names_from = "sero", values_from = "IFR_cols") %>%
  dplyr::ungroup(.) %>%
  dplyr::rename(seroprev_reg = reg,
                seroprev_serorev = serorev)



#......................
# get crude IFRs
#......................
set.seed(48)
source("R/monte_carlo_cis.R")

# calculate CIs for binomial
crude_IFRs_binomial <- dsc_agedat %>%
  dplyr::select(c("study_id", "plotdat")) %>%
  dplyr::filter(!study_id %in% c("ITA1", "SWE1", "DNK1")) %>%
  tidyr::unnest(cols = plotdat) %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::filter(obsday == seromidpt) %>% # latest serostudy
  dplyr::ungroup(.) %>%
  dplyr::group_by(study_id, ageband) %>%
  dplyr::summarise(cumdeaths = sum(cumdeaths),
                   popn = sum(popn),
                   n_positive = sum(n_positive),
                   n_tested = sum(n_tested)) %>%
  dplyr::select(c("study_id", "ageband", "cumdeaths", "popn", "n_positive", "n_tested")) %>%
  dplyr::group_by(study_id, ageband) %>%
  dplyr::mutate(seroprev = n_positive/n_tested,
                ifr_range = purrr::map(cumdeaths, get_binomial_monte_carlo_cis, popN = popn,
                                       npos = n_positive, ntest = n_tested, iters = 1e5),
                crudeIFR = cumdeaths/((seroprev * popn) + cumdeaths),
                lower_ci = purrr::map_dbl(ifr_range, quantile, 0.025),
                upper_ci = purrr::map_dbl(ifr_range, quantile, 0.975)) %>%
  dplyr::select(c("study_id", "ageband", "crudeIFR", "lower_ci", "upper_ci")) %>%
  dplyr::ungroup(.)

# calculate CIs for logit
crude_IFRs_logit <- dsc_agedat %>%
  dplyr::select(c("study_id", "plotdat")) %>%
  dplyr::filter(study_id %in% c("ITA1", "SWE1", "DNK1")) %>%
  tidyr::unnest(cols = plotdat) %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::filter(obsday == seromidpt) %>% # latest serostudy
  dplyr::group_by(study_id, ageband) %>%
  dplyr::summarise(cumdeaths = sum(cumdeaths),
                   popn = sum(popn),
                   seroprev = mean(seroprev),
                   serolci = mean(serolci),
                   serouci = mean(serouci)) %>%
  dplyr::select(c("study_id", "ageband", "cumdeaths", "popn", "seroprev",  "serolci", "serouci")) %>%
  dplyr::group_by(study_id, ageband) %>%
  dplyr::mutate(SE = (COVIDCurve:::logit(serouci) - COVIDCurve:::logit(serolci))/(1.96 * 2))  %>%
  dplyr::mutate(ifr_range = purrr::map(cumdeaths, get_normal_monte_carlo_cis, popN = popn,
                                       mu = seroprev, sigma = SE, iters = 1e5),
                crudeIFR = cumdeaths/((seroprev * popn) + cumdeaths),
                lower_ci = purrr::map_dbl(ifr_range, quantile, 0.025),
                upper_ci = purrr::map_dbl(ifr_range, quantile, 0.975)) %>%
  dplyr::select(c("study_id", "ageband", "crudeIFR", "lower_ci", "upper_ci"))

# out
crudeIFRs_CI <- dplyr::bind_rows(crude_IFRs_binomial, crude_IFRs_logit)
crude_IFRs_column <- crudeIFRs_CI %>%
  dplyr::mutate(crudeIFR = round(crudeIFR * 100, 2),
                lower_ci = round(lower_ci * 100, 2),
                upper_ci = round(upper_ci * 100, 2),
                crude_IFRs = paste0(crudeIFR, " (", lower_ci, ", ", upper_ci, ")"))

#......................
# make final table
#......................
dplyr::left_join(simp_seroprevdat, seroprev_age_column, by = c("study_id", "ageband")) %>%
  dplyr::left_join(., crude_IFRs_column, by = c("study_id", "ageband")) %>%
  dplyr::left_join(., modIFR_age_column, by = c("study_id", "ageband")) %>%
  dplyr::select(c("location", "study_id", "sero_midday", "obs_input_seropev", "ageband", "seroprev_reg", "seroprev_serorev", "crude_IFRs", "regIFR", "serorevIFR")) %>%
  dplyr::left_join(., order, by = "study_id") %>%
  dplyr::mutate(age_low = as.numeric(stringr::str_extract(ageband, "[0-9]+(?=\\,)"))) %>%
  dplyr::arrange(order, age_low) %>%
  dplyr::select(-c("order", "age_low")) %>%
  dplyr::select("location", "ageband", "sero_midday", "obs_input_seropev", "seroprev_reg", "seroprev_serorev", "crude_IFRs", "regIFR", "serorevIFR") %>%  # fix order
  readr::write_tsv(., path = "tables/final_tables/age_specific_ifr_data.tsv")

#............................................................
#---- Best Age IFR Est Table  #----
#...........................................................
#......................
# unpack data for est
#......................
inputdat <- retmapIFR %>%
  dplyr::select(c("study_id", "sero", "strataIFRret")) %>%
  tidyr::unnest(cols = "strataIFRret") %>%
  dplyr::rename(ageband = strata) %>%
  dplyr::mutate(age_mid = purrr::map_dbl(ageband, get_mid_age)) %>%
  dplyr::select(c("study_id", "sero", "ageband", "age_mid", "mean")) %>%
  dplyr::rename(age = age_mid,
                mu = mean) %>%
  dplyr::group_by(sero) %>%
  tidyr::nest(.)

#......................
# get log variance for each age band
#......................
get_age_log_variance <- function(path) {
  # read in
  IFRmodel_inf <- readRDS(path)
  # get variance
  var <- IFRmodel_inf$mcmcout$output %>%
    dplyr::filter(stage == "sampling" & rung == "rung1") %>%
    dplyr::select(c("iteration", dplyr::starts_with("ma"))) %>%
    tidyr::pivot_longer(., cols = -c("iteration"), names_to = "param", values_to = "est") %>%
    dplyr::group_by(param) %>%
    dplyr::mutate(est = log(est)) %>%
    dplyr::summarise(param_logvar = var(est))
  # get data dictionary
  dict <- get_data_dict(path) %>%
    dplyr::rename(param = Strata)
  # out
  dplyr::left_join(dict, var, by = "param") %>%
    dplyr::select(-c("param"))
}
logvardat <- retmapIFR %>%
  dplyr::select(c("study_id", "path", "sero")) %>%
  dplyr::mutate(paramvar = purrr::map(path, get_age_log_variance)) %>%
  dplyr::select(-c("path")) %>%
  tidyr::unnest(cols = "paramvar") %>%
  dplyr::group_by(sero) %>% # put back in structure to match ifr data
  tidyr::nest(.) %>%
  dplyr::rename(vardata = data)


#..........................................
# Fitting models to mu and residuals
# for no seroreversion
#..........................................
fit_pred_intervals_log_linear_mod <- function(ifrdata, vardata) {
  # join age-specific ifrs across studies w/ variance for each age group
  moddat <- dplyr::left_join(ifrdata, vardata, by = c("study_id", "ageband"))
  # fit a weighted log-linear model
  logmod <- lm(log(mu) ~ age, data = moddat, weights = 1/param_logvar)

  # get residual variance in 10-year age groups
  moddat$residuals <- logmod$residuals
  logmodcoeffs <- logmod$coefficients

  # internal check for leverage
  # plot(moddat$age, moddat$residuals)
  moddat$cut <- cut(moddat$age, breaks = c(0, seq(9, 89, 10), 999))
  new_dat <- data.frame(ageband = levels(moddat$cut))
  new_dat$var <- mapply(var, split(moddat$residuals, f = moddat$cut))

  # get first- and second-order age terms
  new_dat$age <- purrr::map_dbl(new_dat$ageband, get_mid_age)
  new_dat$age_squared <- new_dat$age^2

  # fit quadratic model to log of residual variance
  mod2 <- lm(log(var) ~ age + age_squared, data = new_dat)
  new_dat$var_fit <- exp(predict(mod2, newdata = new_dat))

  # get corresponding prediction from optim fit
  new_dat$fit <- logmodcoeffs[2] * new_dat$age + logmodcoeffs[1]

  # get prediction intervals in linear space
  new_dat_log <- new_dat
  new_dat_log$Q025 <- qnorm(0.025, mean = new_dat$fit, sd = sqrt(new_dat$var_fit))
  new_dat_log$Q20 <- qnorm(0.2, mean = new_dat$fit, sd = sqrt(new_dat$var_fit))
  new_dat_log$Q50 <- qnorm(0.5, mean = new_dat$fit, sd = sqrt(new_dat$var_fit))
  new_dat_log$Q80 <- qnorm(0.8, mean = new_dat$fit, sd = sqrt(new_dat$var_fit))
  new_dat_log$Q975 <- qnorm(0.975, mean = new_dat$fit, sd = sqrt(new_dat$var_fit))

  # get prediction intervals in log space
  new_dat_linear <- new_dat
  new_dat_linear$Q025 <- qlnorm(0.025, mean = new_dat$fit, sd = sqrt(new_dat$var_fit))
  new_dat_linear$Q20 <- qlnorm(0.2, mean = new_dat$fit, sd = sqrt(new_dat$var_fit))
  new_dat_linear$mean = exp(new_dat$fit + new_dat$var_fit/2)
  new_dat_linear$Q50 <- qlnorm(0.5, mean = new_dat$fit, sd = sqrt(new_dat$var_fit))
  new_dat_linear$Q80 <- qlnorm(0.8, mean = new_dat$fit, sd = sqrt(new_dat$var_fit))
  new_dat_linear$Q975 <- qlnorm(0.975, mean = new_dat$fit, sd = sqrt(new_dat$var_fit))

  # out
  out <- list(residuals = moddat[, c("age", "residuals")],
              newvar = new_dat[, c("age", "var")],
              linear_new_dat = new_dat_linear,
              log_new_data = new_dat_log,
              logmod = logmod,
              varmod = mod2)
  return(out)
}

#.............
# get fits
#.............
ifrdat <- inputdat %>%
  rename(ifrdata = data) %>%
  dplyr::left_join(., logvardat) %>%
  dplyr::mutate(bestestmod = purrr::map2(ifrdata, vardata, fit_pred_intervals_log_linear_mod),
                logmod = purrr::map(bestestmod, "logmod"),
                varmod = purrr::map(bestestmod, "varmod"),
                residuals = purrr::map(bestestmod, "residuals"),
                logvar = purrr::map(bestestmod, "newvar"),
                linear_predints = purrr::map(bestestmod, "linear_new_dat"),
                log_predints = purrr::map(bestestmod, "log_new_data"))

ifrdat$logmod[[1]]
ifrdat$logmod[[2]]


#......................
# internal checks
#......................
ggplot() +
  geom_ribbon(data = ifrdat$log_predints[[1]], aes(x = age, ymin = Q025, ymax = Q975), alpha = 0.75) +
  geom_point(data = ifrdat$ifrdata[[1]], aes(x = age, y = log(mu))) +
  xlim(c(0,100))

ggplot() +
  geom_ribbon(data = ifrdat$log_predints[[2]], aes(x = age, ymin = Q025, ymax = Q975), alpha = 0.75) +
  geom_point(data = ifrdat$ifrdata[[2]], aes(x = age, y = log(mu))) +
  xlim(c(0,100))


#......................
# plot out residuals
#......................
pA <- ifrdat$residuals[[1]] %>%
  ggplot() +
  geom_point(aes(x = age, y = residuals)) +
  xlab("Age") + ylab("Residuals") +
  xyaxis_plot_theme

pB <- ifrdat$residuals[[2]] %>%
  ggplot() +
  geom_point(aes(x = age, y = residuals)) +
  xlab("Age") + ylab("Residuals") +
  xyaxis_plot_theme

jpeg("figures/final_figures/bestest_residuals_overall_ifr_est.jpg",
     width = 8, heigh = 6, units = "in", res = 500)
cowplot::plot_grid(pA, pB, nrow = 1, labels = c("(A)", "(B)"))
graphics.off()

#......................
# send out variance of new dat for supp
#......................
t1 <- ifrdat$logvar[[1]]
t2 <- ifrdat$logvar[[2]]
dplyr::left_join(t1, t2, by = "age") %>%
  magrittr::set_colnames(c("age", "noserorev_var", "serorev_var")) %>%
  dplyr::mutate_if(is.numeric, round, 4) %>%
  readr::write_tsv(., "tables/final_tables/bestest_ifr_variance_for_mod.tsv")


#......................
# write out the age-specific best fit
#......................
ifrdat %>%
  tidyr::unnest(cols = "linear_predints") %>%
  dplyr::select(c("sero", "ageband", "Q025", "Q50", "Q975")) %>%
  dplyr::mutate(Q025 = round(Q025 * 100, 2),
                Q50 = round(Q50 * 100, 2),
                Q975 = round(Q975 * 100, 2)) %>%
  dplyr::mutate(bestest = paste0(Q50, " (", Q025, ", ", Q975, ")")) %>%
  dplyr::select(c("sero", "ageband", "bestest")) %>%
  tidyr::pivot_wider(., names_from = "sero", values_from = "bestest") %>%
  readr::write_tsv(., path = "tables/final_tables/overall_best_est_for_age_IFRs.tsv")
# look at mean for reference
ifrdat %>%
  tidyr::unnest(cols = "linear_predints") %>%
  dplyr::select(c("sero", "ageband", "mean"))

#............................................................
#---- Best Overall IFR Est Rows  #----
#...........................................................
#......................
# tidy up UN World Population Prospects
#........................
wpp <- readxl::read_excel("data/raw/wpp_un_agepopulations.xlsx")
wpp <- wpp[wpp$`Reference date (as of 1 July)` == 2020,]
wpp <- wpp %>%
  dplyr::select(c("Region, subregion, country or area *", "0-4", "5-9", "10-14", "15-19", "20-24",
                  "25-29", "30-34", "35-39", "40-44", "50-54", "55-59", "60-64", "65-69", "70-74",
                  "75-79", "80-84", "85-89", "90-94", "95-99", "100+"))
colnames(wpp)[1] <- "georegion"

wpp <- wpp %>%
  dplyr::filter(georegion %in% c("Madagascar", "Nicaragua", "Grenada", "Malta")) %>%
  tidyr::pivot_longer(., cols = -c("georegion"), names_to = "ageband", values_to = "popN") %>%
  dplyr::mutate(popN = as.numeric(popN)*1e3) # wpp adjustment

#......................
# cuts
#......................
wpp <- wpp %>%
  dplyr::mutate(
    ageband = ifelse(ageband == "100+", "95-99", ageband),
    age_high = as.numeric(stringr::str_split_fixed(ageband, "-", n=2)[,2]),
    ageband = cut(age_high,
                  breaks = c(0, seq(9, 89, by = 10), 999)),
    ageband = as.character(ageband)) %>%
  dplyr::group_by(georegion, ageband) %>%
  dplyr::summarise(pop_size = sum(popN)) %>%
  dplyr::mutate(age_high = as.numeric(stringr::str_extract(ageband, "[0-9]+?(?=])"))) %>%
  dplyr::arrange(age_high) %>%
  dplyr::select(-c("age_high")) %>%
  dplyr::group_by(georegion) %>%
  tidyr::nest(.) %>%
  dplyr::rename(popN = data)


#........................
# Monte Carlo Calculations
#........................
calc_monte_carlo_overall_IFRs <- function(new_dat, popN, reps = 1e5) {
  if (!all(c("ageband", "pop_size") %in% colnames(popN))) {
    stop()
  }
  # bring in demog
  new_dat <- dplyr::left_join(new_dat, popN)

  # calculate
  z <- mapply(function(i) {
    sum(new_dat$pop_size * rlnorm(nrow(new_dat), meanlog = new_dat$fit, sdlog = sqrt(new_dat$var_fit)))
  }, seq_len(reps))
  z <- z / sum(new_dat$pop_size)

  # get prediction intervals
  Q025 <- round(quantile(z, 0.025) * 100, 2)
  Q50 <- round(quantile(z, 0.50) * 100, 2)
  Q975 <- round(quantile(z, 0.975) * 100, 2)
  # get mean
  mean <- exp(new_dat$fit + new_dat$var_fit/2)
  mean <- sum( mean * (new_dat$pop_size/sum(new_dat$pop_size)) )
  mean <- round(mean * 100, 2)

  # out
  ret <- tibble::tibble(Q025 = Q025,
                        Q50 = Q50,
                        mean = mean,
                        Q975 = Q975)
  return(ret)
}

# tidy up pieces
overall_IFR_best_est <- tidyr::expand_grid(ifrdat, wpp) %>%
  dplyr::select(c("sero", "linear_predints", "georegion", "popN")) %>%
  dplyr::rename(new_dat = linear_predints)
overall_IFR_best_est$bestest <- purrr::pmap(overall_IFR_best_est[, c("new_dat", "popN")],
                                            calc_monte_carlo_overall_IFRs,
                                            reps = 1e5)

# send out
overall_IFR_best_est %>%
  dplyr::select(c("sero", "georegion", "bestest")) %>%
  tidyr::unnest(cols = "bestest") %>%
  dplyr::mutate(bestest = paste0(Q50, " (", Q025, ", ", Q975, ")")) %>%
  dplyr::select(c("sero", "georegion", "bestest")) %>%
  tidyr::pivot_wider(., names_from = "sero", values_from = "bestest") %>%
  readr::write_tsv(., path = "tables/final_tables/overall_best_est_for_georegions_IFRs.tsv")

# look at mean for reference
overall_IFR_best_est %>%
  dplyr::select(c("sero", "georegion", "bestest")) %>%
  tidyr::unnest(cols = "bestest")

#............................................................
#---- Figure of Age Specific Results #----
#...........................................................
#......................
# tidy modelled and plot
#......................
modIFR_age <- modIFR_age %>%
  dplyr::mutate(age_mid = purrr::map_dbl(ageband, get_mid_age))

#......................
# get log transformed variables
#......................
log_retmapIFR <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(strataIFRret_log = purrr::map(path, get_log10_transformed_IFR_cred_intervals, by_chain = F))

log_retmapIFR_dat <- log_retmapIFR %>%
  tidyr::unnest(cols = "strataIFRret_log") %>%
  dplyr::rename(Strata = param) %>%
  dplyr::left_join(., datdict) %>%
  dplyr::mutate(age_mid = purrr::map_dbl(ageband, get_mid_age))

#......................
# Panel A -- No Seroreversion Linear Space
#......................
ribbondat <- ifrdat$linear_predints[ifrdat$sero == "reg"][[1]] %>%
  dplyr::mutate(Q025 = Q025 *100,
                Q20 = Q20 * 100,
                Q80 = Q80 * 100,
                Q975 = Q975 * 100)


PanelA <- modIFR_age %>%
  dplyr::left_join(., locatkey, by = "study_id") %>%
  dplyr::filter(sero == "reg") %>%
  ggplot() +
  geom_ribbon(data = ribbondat,
              aes(x = age, ymin = Q025, ymax = Q975),
              fill = "#d9d9d9", alpha = 0.8) +
  geom_ribbon(data = ribbondat,
              aes(x = age, ymin = Q20, ymax = Q80),
              fill = "#969696", alpha = 0.8) +
  geom_pointrange(aes(x = age_mid, y = median, ymin = LCI, ymax = UCI,
                      color =  location),
                  alpha = 0.75, shape = 16, size = 0.9) +
  scale_color_manual("Location", values = mycolors) +
  ylab("IFR (95% CrI)") + xlab("Age (years)") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(0.05, 0.25, 0.25, 1),"cm"))




#......................
# Panel B --  No Seroreversion Log10 Space
#......................
ribbondat <- ifrdat$log_predints[ifrdat$sero == "reg"][[1]] %>%
  dplyr::mutate(Q025 = Q025/log(10), # tranformation
                Q20 = Q20/log(10),
                Q80 = Q80/log(10),
                Q975 = Q975 /log(10))

PanelB <- log_retmapIFR_dat %>%
  dplyr::left_join(., locatkey, by = "study_id") %>%
  dplyr::filter(sero == "reg") %>%
  ggplot() +
  geom_ribbon(data = ribbondat,
              aes(x = age, ymin = Q025, ymax = Q975),
              fill = "#d9d9d9", alpha = 0.8) +
  geom_ribbon(data = ribbondat,
              aes(x = age, ymin = Q20, ymax = Q80),
              fill = "#969696", alpha = 0.8) +
  geom_pointrange(aes(x = age_mid, y = median, ymin = LCI, ymax = UCI, color =  location),
                  alpha = 0.75, shape = 16, size = 0.9) +
  scale_color_manual("Location", values = mycolors) +
  ylab("Log-10 IFR (95% CrI)") + xlab("Age (years)") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(0.05, 0.25, 0.25, 1),"cm"))

#......................
# Panel C -- Seroreversion Linear Space
#......................
ribbondat <- ifrdat$linear_predints[ifrdat$sero == "serorev"][[1]] %>%
  dplyr::mutate(Q025 = Q025 *100,
                Q20 = Q20 * 100,
                Q80 = Q80 * 100,
                Q975 = Q975 * 100)


PanelC <- modIFR_age %>%
  dplyr::left_join(., locatkey, by = "study_id") %>%
  dplyr::filter(sero == "serorev") %>%
  ggplot() +
  geom_ribbon(data = ribbondat,
              aes(x = age, ymin = Q025, ymax = Q975),
              fill = "#d9d9d9", alpha = 0.8) +
  geom_ribbon(data = ribbondat,
              aes(x = age, ymin = Q20, ymax = Q80),
              fill = "#969696", alpha = 0.8) +
  geom_pointrange(aes(x = age_mid, y = median, ymin = LCI, ymax = UCI,
                      color =  location),
                  alpha = 0.75, shape = 16, size = 0.9) +
  scale_color_manual("Location", values = mycolors) +
  ylab("IFR (95% CrI)") + xlab("Age (years)") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(1, 0.25, 0.25, 1),"cm"))

#......................
# Panel D -- Seroreversion Log10 Space
#......................
ribbondat <- ifrdat$log_predints[ifrdat$sero == "serorev"][[1]] %>%
  dplyr::mutate(Q025 = Q025/log(10), # tranformation
                Q20 = Q20/log(10),
                Q80 = Q80/log(10),
                Q975 = Q975 /log(10))

PanelD <- log_retmapIFR_dat %>%
  dplyr::left_join(., locatkey, by = "study_id") %>%
  dplyr::filter(sero == "serorev") %>%
  ggplot() +
  geom_ribbon(data = ribbondat,
              aes(x = age, ymin = Q025, ymax = Q975),
              fill = "#d9d9d9", alpha = 0.8) +
  geom_ribbon(data = ribbondat,
              aes(x = age, ymin = Q20, ymax = Q80),
              fill = "#969696", alpha = 0.8) +
  geom_pointrange(aes(x = age_mid, y = median, ymin = LCI, ymax = UCI, color =  location),
                  alpha = 0.75, shape = 16, size = 0.9) +
  scale_color_manual("Location", values = mycolors) +
  ylab("Log-10 IFR (95% CrI)") + xlab("Age (years)") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(1, 0.25, 0.25, 1),"cm"))



#......................
# bring together
#......................
PanelAnl <- PanelA + theme(legend.position = "none")
PanelBnl <- PanelB + theme(legend.position = "none")
PanelCnl <- PanelC + theme(legend.position = "none")
PanelDnl <- PanelD + theme(legend.position = "none")
mainFig <- cowplot::plot_grid(PanelAnl, PanelBnl, PanelCnl, PanelDnl,
                              align = "h", ncol = 2, nrow = 2,
                              labels = c("(A)", "", "(B)", ""))
legend <- cowplot::get_legend(PanelA +
                                guides(color = guide_legend(nrow = 3)) +
                                theme(legend.position = "bottom"))
mainFig <- cowplot::plot_grid(mainFig, legend,
                              nrow = 2, ncol = 1, rel_heights = c(0.9, 0.15))

jpeg("figures/final_figures/IFR_age_spec_logplot.jpg",
     width = 11, height = 8, units = "in", res = 600)
plot(mainFig)
graphics.off()




