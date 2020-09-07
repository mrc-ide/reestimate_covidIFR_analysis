#....................................................................................................
## Purpose:
##
## Notes:
#....................................................................................................
library(tidyverse)
library(brms)
source("R/my_themes.R")


#......................
# read data
#......................
serotime <- readr::read_csv("data/raw/convalescent_plasma_longitudinal.csv")
serotime <- serotime %>%
  dplyr::mutate(sex = factor(sex, levels = c("F", "M")),
                days_post_symptoms = as.numeric( lubridate::mdy(date_sample_collected) - lubridate::mdy(date_symptom_onset)),
                months_post_symptoms = days_post_symptoms/30)  # re-center days to something more central

#............................................................
# Explore Data and Summarize Dist
#...........................................................
summary(serotime)

serotime %>%
  ggplot() +
  geom_line(aes(x = days_post_symptoms, y = abbott_sc, group = donor_id)) +
  xyaxis_plot_theme


# sex pattern not apparent
serotime %>%
  ggplot() +
  geom_line(aes(x = days_post_symptoms, y = abbott_sc, group = donor_id)) +
  facet_wrap(. ~ sex) +
  xyaxis_plot_theme


# lose more in hospitalized
serotime %>%
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(diffAb = max(abbott_sc) - min(abbott_sc)) %>%
  ggplot() +
  geom_histogram(aes(x = abbott_sc, y = ..density..)) +
  facet_wrap(. ~ hospitalized) +
  xyaxis_plot_theme


# lose more in hospitalized
serotime %>%
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(diffAb = max(abbott_sc) - min(abbott_sc),
                difftime = max(days_post_symptoms),
                age = mean(age)) %>%
  ggplot() +
  geom_point(aes(x = difftime, y = diffAb, color = age)) +
  scale_color_gradientn("Age",
                        colors = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))) +
  ylab("Diff") +
  facet_wrap(. ~ hospitalized) +
  xyaxis_plot_theme


# function of start
serotime %>%
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(diffAb = log2(min(abbott_sc)/max(abbott_sc)),
                difftime = max(days_post_symptoms),
                age = mean(age)) %>%
  ggplot() +
  geom_point(aes(x = difftime, y = diffAb, color = age)) +
  ylab("Fold Change") +
  scale_color_gradientn("Age",
                        colors = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))) +
  facet_wrap(. ~ hospitalized) +
  xyaxis_plot_theme



#............................................................
# MLModeling
#...........................................................
# lognormal, just intercept
fit.ln.1 <- brm(data = serotime, family = lognormal,
                abbott_sc ~ 1   + (1 | donor_id),
                prior = c(prior(normal(0, 5), class = "Intercept"),
                          prior(cauchy(0, 1), class = "sd")),
                iter = 3000, warmup = 1000, chains = 3, cores = 1,
                seed = 48)
summary(fit.ln.1, waic = TRUE)
# launch_shinystan(fit.ln.1)


# lognormal, just time random intercept
fit.ln.2 <- brm(data = serotime, family = lognormal,
                abbott_sc ~ months_post_symptoms   + (1 | donor_id),
                prior = c(prior(normal(0, 5), class = "Intercept"),
                          prior(normal(0, 5), class = "b", coef = "months_post_symptoms"),
                          prior(cauchy(0, 1), class = "sd")),
                iter = 3000, warmup = 1000, chains = 3, cores = 1,
                seed = 48)
summary(fit.ln.2, waic = TRUE)
# launch_shinystan(fit.ln.2)

# lognormal, time and hospital random intercept
fit.ln.3 <- brm(data = serotime, family = lognormal,
                abbott_sc ~ months_post_symptoms + hospitalized + (1 | donor_id),
                prior = c(prior(normal(0, 5), class = "Intercept"),
                          prior(normal(0, 5), class = "b"),
                          prior(cauchy(0, 1), class = "sd")),
                iter = 3000, warmup = 1000, chains = 3, cores = 1,
                seed = 48)
summary(fit.ln.3, waic = TRUE)
# launch_shinystan(fit.ln.3)

# time random slope
fit.ln.4 <- brm(data = serotime, family = lognormal,
                abbott_sc ~ 1 + (months_post_symptoms | donor_id),
                prior = c(prior(normal(0, 5), class = "Intercept"),
                          prior(cauchy(0, 1), class = "sd")),
                iter = 3000, warmup = 1000, chains = 3, cores = 1,
                seed = 48)
summary(fit.ln.4, waic = TRUE)
# launch_shinystan(fit.ln.4)


# time random slope, fixed hospital
fit.ln.5 <- brm(data = serotime, family = lognormal,
                abbott_sc ~ hospitalized + (months_post_symptoms | donor_id),
                prior = c(prior(normal(0, 5), class = "Intercept"),
                          prior(normal(0, 5), class = "b"),
                          prior(cauchy(0, 1), class = "sd")),
                iter = 2000, warmup = 1000, chains = 3, cores = 1,
                seed = 48)
summary(fit.ln.5, waic = TRUE)
# launch_shinystan(fit.ln.5)


# compare
fit.ln.1 <- add_criterion(fit.ln.1, "waic")
fit.ln.1 <- add_criterion(fit.ln.1, "loo")
fit.ln.2 <- add_criterion(fit.ln.2, "waic")
fit.ln.2 <- add_criterion(fit.ln.2, "loo")
fit.ln.3 <- add_criterion(fit.ln.3, "waic")
fit.ln.3 <- add_criterion(fit.ln.3, "loo")
fit.ln.4 <- add_criterion(fit.ln.4, "waic")
fit.ln.4 <- add_criterion(fit.ln.4, "loo")
fit.ln.5 <- add_criterion(fit.ln.5, "waic")
fit.ln.5 <- add_criterion(fit.ln.5, "loo")




loo_compare(fit.ln.1, fit.ln.2, fit.ln.3, fit.ln.4, fit.ln.5, criterion = "loo") %>%
  print(simplify = F)
loo_compare(fit.ln.1, fit.ln.2, fit.ln.3, fit.ln.4, fit.ln.5, criterion = "waic") %>%
  print(simplify = F)

model_weights(fit.ln.1, fit.ln.2, fit.ln.3, fit.ln.4, fit.ln.5,
              weights = "waic") %>%
  round(digits = 2)


#............................................................
# Posterior Predictions from best model
#...........................................................
donorids <- unique(serotime$donor_id)
newdat <- lapply(donorids, function(x){
  tibble::tibble(donor_id = x,
                 months_post_symptoms = seq(min(serotime$days_post_symptoms)/30, max(serotime$days_post_symptoms)/30,
                                            length.out = 10))}) %>%
  dplyr::bind_rows()

preds <- fitted(fit.ln.4, newdata = newdat, type = "response")
preds <- cbind.data.frame(newdat, preds)
preds %>%
  ggplot() +
  geom_line(aes(x = months_post_symptoms, y = Estimate, group = factor(donor_id)), alpha = 0.8, color = "#bdbdbd") +
  geom_line(data = serotime, aes(x = months_post_symptoms, y = abbott_sc, group = donor_id), color = "#0000F9", alpha = 0.5) +
  xyaxis_plot_theme +
  ylab("Pred Abott SC") +
  xlab("Months Post-Sxs")




#............................................................
# optimization function to find posterior point
# that is closest to the time of cutoff for seroreversion
# based on a prior assumption of best model fit
#...........................................................

# internal function
find_cutoff_response_time_per_donor <- function(cutoff, donorrow, mod) {
  if (!all(colnames(donorrow) %in% c("donor_id", "months_post_symptoms")) ) {
    stop("Assumed we were fitting the model with only the donor_id and months_post_sx variables.
         This function is not generalizable")
  }
  # internal cost function
  cost <- function(par, mod, dat, cutoff) {
    dat[,"months_post_symptoms"] <- par
    pred <- fitted(mod, newdata = dat, type = "response")[1, "Estimate"]
    cost <- (cutoff - pred)^2
    return(cost)
  }
  out <- optim(par = 1, fn = cost, mod = mod, dat = donorrow, cutoff = cutoff,
               method = "L-BFGS-B", lower = .Machine$double.xmin)$par
  return(out)
}

# split data for run
donor_id_dat <- tibble::tibble(
  donor_id = unique(serotime$donor_id),
  months_post_symptoms = 0)
donor_id_dat <- split(donor_id_dat, 1:nrow(donor_id_dat))

# run
serorevert_times <- furrr::future_map_dbl(donor_id_dat, find_cutoff_response_time_per_donor,
                                          cutoff = 1.4, mod = fit.ln.2)
serorevert_times <- serorevert_times * 30 # back to days
#............................................................
# fitdistr
#...........................................................
wb_fit <- fitdistrplus::fitdist(serorevert_times, distr = "weibull")
gm_fit <- fitdistrplus::fitdist(serorevert_times, distr = "gamma")
lg_fit <- fitdistrplus::fitdist(serorevert_times, distr = "lnorm")
exp_fit <- fitdistrplus::fitdist(serorevert_times, distr = "exp")

# compare
par(mfrow = c(2, 2))
plot.legend <- c("Weibull", "gamma", "lognormal", "exponential")
fitdistrplus::denscomp(list(wb_fit, gm_fit, lg_fit, exp_fit), legendtext = plot.legend)
fitdistrplus::qqcomp(list(wb_fit, gm_fit, lg_fit, exp_fit), legendtext = plot.legend)
fitdistrplus::cdfcomp(list(wb_fit, gm_fit, lg_fit, exp_fit), legendtext = plot.legend)
fitdistrplus::ppcomp(list(wb_fit, gm_fit, lg_fit, exp_fit), legendtext = plot.legend)
# gof stat
fitdistrplus::gofstat(list(wb_fit, gm_fit, lg_fit, exp_fit),
        fitnames =  c("Weibull", "gamma", "lognormal", "exponential"))
# weibull best
wb_fit_boot <- fitdistrplus::bootdist(wb_fit, niter = 1e3)
summary(wb_fit_boot)
plot(wb_fit_boot)
par(mfrow = c(1, 2))
plot(density(wb_fit_boot$estim$shape))
plot(density(wb_fit_boot$estim$scale))
# looks like normal distributions to me
# just for fun what does the sd look like
fitdistrplus::fitdist(wb_fit_boot$estim$scale, "norm")
fitdistrplus::fitdist(wb_fit_boot$estim$shape, "norm")
