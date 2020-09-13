#....................................................................................................
## Purpose: Determine the distribution for the onset-of-symptoms to seroreversion
##
## Notes: Data shared from Muecksch et. al 2020
#....................................................................................................
library(tidyverse)
library(brms)
library(survival)
library(survminer)
source("R/my_themes.R")
source("R/extra_plotting_functions.R")
set.seed(48)

#......................
# read data
#......................
serotime <- readxl::read_excel("data/raw/shared/2020_08_06_SR1407_meta_analysis.xlsx", sheet = 2) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  magrittr::set_colnames(gsub(" ", "_", colnames(.))) %>%
  magrittr::set_colnames(gsub("/", "_", colnames(.))) %>%
  magrittr::set_colnames(gsub("≥1", "gt1", colnames(.))) %>%
  dplyr::rename(sex = gender,
                donor_id = sr1407,
                days_post_symptoms = post_sx_days,
                roche_gt1 = `roche_≥1`) %>%
  dplyr::mutate(sex = factor(sex, levels = c("F", "M")),
                donor_id = factor(donor_id),
                months_post_symptoms = days_post_symptoms/30)    # re-center days to months
# take to long format
serotime <- serotime %>%
  dplyr::select(-c("siemens")) %>%
  tidyr::pivot_longer(., cols = c("abbott_s_c", "diasorin_au_ml", "siemens_no.", "roche_gt1"),
                      names_to = "assay", values_to = "titres") %>%
  dplyr::mutate(assay = stringr::str_split_fixed(assay, "_", n = 2)[,1],
                assay = factor(assay, levels = c("diasorin", "siemens", "abbott", "roche"),
                               labels = c("Diasorin", "Siemens", "Abbott", "Roche")))

# drop hospitalized
serotime <- serotime %>%
  dplyr::filter(hosp == "N")

#............................................................
#---- Explore Data and Summarize Dist #----
#...........................................................
summary(serotime)
# must be positive at baseline
thresholds <- tibble::tibble(assay = c("Diasorin", "Siemens", "Abbott", "Roche"),
                             threshold = c(15, 1, 1.4, 1))
# look at post pcr -- figure 1
serotime %>%
  dplyr::left_join(., thresholds, by = "assay") %>%
  dplyr::filter(assay == "Abbott") %>%
  ggplot() +
  geom_line(aes(x = post_pcr_days, y = titres, group = donor_id, color = assay)) +
  geom_hline(aes(yintercept = threshold), linetype = "dashed") +
  scale_color_manual(values = c("#3182bd", "#31a354", "#de2d26", "#756bb1")) +
  facet_wrap(~assay, scales = "free_y") +
  xyaxis_plot_theme
# look at post sx instead of pcr (but now keep everyone)
serotime %>%
  ggplot() +
  geom_line(aes(x = days_post_symptoms, y = titres, group = donor_id, color = assay)) +
  scale_color_manual(values = c("#3182bd", "#31a354", "#de2d26", "#756bb1")) +
  facet_wrap(~assay, scales = "free_y") +
  ylab("Ab. Titres") + xlab("Days Post-Symptom Onset") +
  xyaxis_plot_theme

# sex pattern not apparent
serotime %>%
  ggplot() +
  geom_line(aes(x = days_post_symptoms, y = titres, group = donor_id, color = assay)) +
  scale_color_manual(values = c("#3182bd", "#31a354", "#de2d26", "#756bb1")) +
  facet_grid(assay ~ sex, scales = "free_y") +
  xyaxis_plot_theme


# very few in hospital but may have higher starting baseline
serotime %>%
  ggplot() +
  geom_line(aes(x = days_post_symptoms, y = titres, group = donor_id, color = assay)) +
  scale_color_manual(values = c("#3182bd", "#31a354", "#de2d26", "#756bb1")) +
  facet_grid(assay ~ hosp, scales = "free_y") +
  xyaxis_plot_theme

# function of start -- most decline
serotime %>%
  dplyr::group_by(donor_id, assay) %>%
  dplyr::mutate(diffAb = log2(min(titres)/max(titres)),
                difftime = max(days_post_symptoms)) %>%
  ggplot() +
  geom_point(aes(x = difftime, y = diffAb, color = assay)) +
  scale_color_manual(values = c("#3182bd", "#31a354", "#de2d26", "#756bb1")) +
  facet_wrap(~assay, scales = "free_y") +
  ylab("Fold Change") +
  xyaxis_plot_theme

# function of start and end important by age?
serotime %>%
  dplyr::group_by(donor_id, assay) %>%
  dplyr::mutate(diffAb = log2(min(titres)/max(titres)),
                difftime = max(days_post_symptoms),
                age = mean(age)) %>%
  ggplot() +
  geom_point(aes(x = difftime, y = diffAb, color = age)) +
  ylab("Fold Change") +
  scale_color_gradientn("Age",
                        colors = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))) +
  facet_wrap(~assay, scales = "free_y") +
  xyaxis_plot_theme

#............................................................
#---- Multilevel Modeling #----
#...........................................................
#......................
# data wrangling for modeling
#......................
# subsetting just to the abbott assay
serotime <- serotime %>%
  dplyr::filter(assay == "Abbott") %>%
  dplyr::filter(!is.na(days_post_symptoms))

# must have at least three timepoints (since can fit any line through two)
# sero_nobs <- serotime %>%
#   dplyr::group_by(donor_id) %>%
#   dplyr::summarise(nobs = sum(!is.na(titres)))
# sero_sub <- dplyr::left_join(serotime, sero_nobs, by = c("donor_id")) %>%
#   dplyr::filter(nobs >= 3) %>%
#   dplyr::select(-c("nobs"))

# must be positive at baseline
thresholds <- tibble::tibble(assay = c("Diasorin", "Siemens", "Abbott", "Roche"),
                             threshold = c(15, 1, 1.4, 1))

neg_at_baseline <-  serotime %>%
  dplyr::group_by(donor_id, assay) %>%
  dplyr::filter(days_post_symptoms == min(days_post_symptoms, na.rm = TRUE)) %>%
  dplyr::left_join(., thresholds, by = c("assay")) %>%
  dplyr::mutate(neg_at_baseline = titres < threshold)

neg_at_baseline %>%
  ungroup(.) %>%
  dplyr::filter(neg_at_baseline == T) %>%
  dplyr::select(c("assay", "donor_id")) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::group_by(assay) %>%
  dplyr::summarise(
    n_donors = dplyr::n()
  )
neg_at_baseline <- neg_at_baseline %>%
  dplyr::ungroup(.) %>%
  dplyr::filter(neg_at_baseline == T) %>%
  dplyr::select(c("assay", "donor_id")) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::mutate(drop = TRUE)

# drop those individuals negative at baseline
sero_sub_final <- serotime %>%
  dplyr::left_join(., neg_at_baseline, by = c("assay", "donor_id")) %>%
  dplyr::filter(is.na(drop)) %>%
  dplyr::select(-c("drop"))


#......................
# look at included
#......................
# look at followup
unlft <- sero_sub_final %>%
  dplyr::select(c("donor_id", "assay")) %>%
  dplyr::filter(!duplicated(.))
xtabs(~assay, data = unlft)

fuptims <- sero_sub_final %>%
  dplyr::group_by(donor_id) %>%
  dplyr::summarise(max_obs_sx = max(days_post_symptoms, na.rm = TRUE),
                   max_obs_pcr = max(post_pcr_days, na.rm = TRUE))
summary(fuptims)

# look at characteristics of included
incld <- sero_sub_final %>%
  dplyr::select(c("donor_id", "sex", "age", "hosp")) %>%
  dplyr::filter(!duplicated(.))
summary(incld)
table(sero_sub_final$sex)
table(sero_sub_final$hosp)


# drop to information that we need
sero_sub_final <- sero_sub_final %>%
  dplyr::select(c("donor_id", "assay", "titres", "days_post_symptoms")) %>%
  dplyr::filter(!is.na(titres))


# quick look at finals
sero_sub_final %>%
  dplyr::left_join(., thresholds, by = "assay") %>%
  ggplot() +
  geom_line(aes(x = days_post_symptoms, y = titres, group = donor_id, color = assay),
            color = "#3182bd") +
  geom_hline(aes(yintercept = threshold), linetype = "dashed") +
  facet_wrap(~assay, scales = "free_y") +
  ylab("Ab. Titres") + xlab("Days Post-Symptom Onset") + ggtitle("Final Included") +
  xyaxis_plot_theme

# save out
dir.create("results/sero_reversion/", recursive = T)
saveRDS(sero_sub_final, file = "results/sero_reversion/sero_reversion_incld_data.RDS")


#......................
# modeling functions
#......................
base_model <- function(dat) {
  brms::brm(data = dat, family = lognormal,
            titres ~ 1,
            prior = c(prior(normal(0, 5), class = "Intercept")),
            iter = 1e4, warmup = 5e3, chains = 5, cores = 1,
            seed = 48)
}

RE_intercept_only_model <- function(dat) {
  brms::brm(data = dat, family = lognormal,
            titres ~ 1   + (1 | donor_id),
            prior = c(prior(normal(0, 5), class = "Intercept"),
                      prior(cauchy(0, 1), class = "sd")),
            iter = 1e4, warmup = 5e3, chains = 5, cores = 1,
            seed = 48)
}

RE_sxtime_model <- function(dat) {
  brms::brm(data = dat, family = lognormal,
            titres ~ months_post_symptoms   + (1 | donor_id),
            prior = c(prior(normal(0, 5), class = "Intercept"),
                      prior(normal(0, 5), class = "b", coef = "months_post_symptoms"),
                      prior(cauchy(0, 1), class = "sd")),
            iter = 1e4, warmup = 5e3, chains = 5, cores = 1,
            seed = 48)
}

#......................
# modeling
#......................
#note sigma here is just for constant
#https://discourse.mc-stan.org/t/understanding-a-simple-brms-model-using-the-make-stancode-function/5505/2

sero_sub_mods <- sero_sub_final %>%
  dplyr::group_by(assay) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    basemod = furrr::future_map(data, base_model),
    interceptonly = furrr::future_map(data, RE_intercept_only_model),
    sxtime = furrr::future_map(data, RE_sxtime_model)
  )


#......................
# model convergence
#......................
sero_sub_mods <- sero_sub_mods %>%
  tidyr::pivot_longer(., cols = -c("assay", "data"),
                      names_to = "modlvl", values_to = "mod") %>%
  dplyr::mutate(rhat = purrr::map(mod, brms::rhat))

checkrhat <- sero_sub_mods %>%
  dplyr::filter(modlvl == "sxtime") %>%
  dplyr::select(c("assay", "rhat")) %>%
  unnest(cols = rhat)

#......................
# model comparison
#......................
sero_sub_mods <- sero_sub_mods %>%
  dplyr::mutate(mod = purrr::map(mod, brms::add_criterion, "waic"),
                mod = purrr::map(mod, brms::add_criterion, "loo"))

# abbott
brms::loo_compare(sero_sub_mods$mod[[1]], sero_sub_mods$mod[[2]],
                  sero_sub_mods$mod[[3]],
                  criterion = "loo") %>%
  print(simplify = F)
brms::model_weights(sero_sub_mods$mod[[1]], sero_sub_mods$mod[[2]],
                    sero_sub_mods$mod[[3]],
                    weights = "waic") %>%
  round(digits = 2)

#  models favor RE w/ symptoms (no surprise)

# save out
saveRDS(sero_sub_mods, file = "results/sero_reversion/sero_reversion_model_fits.RDS")

#......................
# model results
#......................
sero_sub_mods_sxtime <- sero_sub_mods %>%
  dplyr::filter(modlvl == "sxtime")
summary(sero_sub_mods_sxtime$mod[[1]])
#............................................................
#---- Extrapolation #----
#...........................................................
donor_ids <- sero_sub_final %>%
  dplyr::group_by(assay) %>%
  tidyr::nest(.) %>%
  dplyr::mutate(donor_ids = map(data, function(x){unique(x$donor_id)})) %>%
  dplyr::ungroup(.) %>%
  dplyr::select(-c("data"))

get_new_dat <- function(data, donor_ids){
  lapply(donor_ids, function(x){ # donor ids from global
    tibble::tibble(donor_id = x,
                   months_post_symptoms = seq(min(data$months_post_symptoms), max(data$months_post_symptoms),
                                              length.out = 10))}) %>%
    dplyr::bind_rows()}

#......................
# Posterior Predictions from best model
#......................
post_mods <- sero_sub_mods %>%
  dplyr::filter(modlvl == "sxtime") %>%
  dplyr::left_join(., donor_ids, by = "assay")
post_mods$newdat <- purrr::pmap(post_mods[, c("donor_ids", "data")], get_new_dat)
# now use model fits
get_preds <- function(mod, newdat) {
  ret <- fitted(mod, newdata = newdat, type = "response")
  ret <- cbind.data.frame(newdat, ret)
  return(ret)
}
post_mods$preds <- purrr::pmap(post_mods[, c("mod", "newdat")], get_preds)

# save out
saveRDS(post_mods, file = "results/sero_reversion/Seroreversion_posterior_interpolations.RDS")



#............................................................
# optimization function to find posterior point
# that is closest to the time of cutoff for seroreversion
# based on a prior assumption of best model fit
#...........................................................
# internal function
find_cutoff_response_time_per_donor <- function(threshold, donor_ids, mod) {
  # Assumed we were fitting the model with only the donor_id and months_post_sx variables
  # This function is not generalizable
  # internal cost function
  cost <- function(par, mod, donor_row, threshold) {
    dat <- tibble::tibble(donor_id = donor_row,
                          months_post_symptoms = par)
    pred <- fitted(mod, newdata = dat, type = "response")[1, "Estimate"]
    cost <- (threshold - pred)^2
    return(cost)
  }
  out <- optim(par = 1, fn = cost, mod = mod, donor_row = donor_ids, threshold = threshold,
               method = "L-BFGS-B", lower = .Machine$double.xmin)$par
  return(out)
}

# get cutoffs
post_mods_optim <- dplyr::left_join(post_mods, thresholds, by = "assay") %>%
  dplyr::select(c("assay", "mod", "donor_ids", "threshold")) %>%
  tidyr::unnest(cols = "donor_ids") %>%
  dplyr::ungroup(.)

post_mods_optim$serorevert_times <- furrr::future_pmap_dbl(post_mods_optim[, c("mod", "donor_ids", "threshold")],
                                                           find_cutoff_response_time_per_donor)

# save out
post_mods_optim %>%
  dplyr::select(c("donor_ids", "assay", "serorevert_times")) %>%
  saveRDS(., file = "results/sero_reversion/sero_reversion_optim_times_extrapolated.RDS")


#............................................................
#---- Fit Distributions #----
#...........................................................
# tidy out
post_mods_optim_fits <- post_mods_optim %>%
  dplyr::mutate(serorevert_times = serorevert_times * 30)  # convert to days for seroreversion framework

summary(post_mods_optim_fits$serorevert_times)




#......................
# fit
#......................
post_mods_optim_fits <- post_mods_optim_fits %>%
  dplyr::group_by(assay) %>%
  tidyr::nest(.) %>%
  dplyr::mutate(wb_fit = purrr::map(data, function(x){fitdistrplus::fitdist(x$serorevert_times, distr = "weibull")}),
                gm_fit = purrr::map(data, function(x){fitdistrplus::fitdist(x$serorevert_times, distr = "gamma")}),
                lg_fit = purrr::map(data, function(x){fitdistrplus::fitdist(x$serorevert_times, distr = "lnorm")}),
                exp_fit = purrr::map(data, function(x){fitdistrplus::fitdist(x$serorevert_times, distr = "exp")}))

# compare
compare_fits <- function(wb_fit, gm_fit, lg_fit, exp_fit) {
  par(mfrow = c(2, 2))
  plot.legend <- c("Weibull", "gamma", "lognormal", "exponential")
  fitdistrplus::denscomp(list(wb_fit, gm_fit, lg_fit, exp_fit), legendtext = plot.legend)
  fitdistrplus::qqcomp(list(wb_fit, gm_fit, lg_fit, exp_fit), legendtext = plot.legend)
  fitdistrplus::cdfcomp(list(wb_fit, gm_fit, lg_fit, exp_fit), legendtext = plot.legend)
  fitdistrplus::ppcomp(list(wb_fit, gm_fit, lg_fit, exp_fit), legendtext = plot.legend)
  # gof stat
  ret <- fitdistrplus::gofstat(list(wb_fit, gm_fit, lg_fit, exp_fit),
                               fitnames =  c("Weibull", "gamma", "lognormal", "exponential"))
  return(ret)
}

post_mods_optim_fits$gof <- purrr::pmap(post_mods_optim_fits[c("wb_fit", "gm_fit", "lg_fit", "exp_fit")],
                                        compare_fits)

# save out
saveRDS(post_mods_optim_fits, file = "results/sero_reversion/sero_reversion_extrapolated_dist_fits.RDS")


#......................
# weibull best
#......................
wb_fit_boot <- fitdistrplus::bootdist(post_mods_optim_fits$wb_fit[[1]], niter = 1e3)
summary(wb_fit_boot)
mixdist::weibullparinv(shape = wb_fit_boot$fitpart$estimate["shape"],
                       scale =  wb_fit_boot$fitpart$estimate["scale"])

plot(wb_fit_boot)
par(mfrow = c(1, 2))
plot(density(wb_fit_boot$estim$shape))
plot(density(wb_fit_boot$estim$scale))
# save out
saveRDS(wb_fit_boot, file = "results/sero_reversion/sero_reversion_abbott_bootstrapped_weibull_params.RDS")

# for priors
summary(wb_fit_boot$estim$shape)
sd(wb_fit_boot$estim$shape)
summary(wb_fit_boot$estim$scale)
sd(wb_fit_boot$estim$scale)



#............................................................
#---- Alternative Survival Analysis #----
#...........................................................
library(survival)
sero_sub_final_survival <- sero_sub_final %>%
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(status = ifelse(titres < 1.4, 1, 0),
                status2=ifelse(min(titres,na.rm=T)<1.4,1,0),
                max_time = ifelse(days_post_symptoms==max(days_post_symptoms),1,0)) %>%
  dplyr::arrange(donor_id,days_post_symptoms)

## those without event. Time1 is final observation time, time2 is missing.
sero_pos<-sero_sub_final_survival %>%
  dplyr::filter(status2==0,days_post_symptoms == max(days_post_symptoms)) %>%
  dplyr::mutate(time1=days_post_symptoms,
                time2=NA)

## those who serorevert
sero_rev<-sero_sub_final_survival %>%
  dplyr::filter(status2==1)
## extract last time observed postiive
sero_rev_time1<- sero_rev %>%
  dplyr::filter(status==0) %>%
  dplyr::filter(days_post_symptoms==max(days_post_symptoms)) %>%
  dplyr::mutate(time1=days_post_symptoms)
#extract first time observed seroreverted (checked that no one becomes positive again after first being negative)
sero_rev_time2<- sero_rev %>%
  dplyr::filter(status==1) %>%
  dplyr::filter(days_post_symptoms==min(days_post_symptoms)) %>%
  dplyr::mutate(time2=days_post_symptoms) %>%
  dplyr::select(donor_id,time2)
## combine serorev
sero_rev_comb<-left_join(sero_rev_time1,sero_rev_time2)

#### combine seropos and serorev
sero_sub_final_survival<-rbind(sero_pos,sero_rev_comb)
sero_sub_final_survival <-sero_sub_final_survival %>%
  dplyr::mutate(time_obs_fail = ifelse(status2==1,time2,time1))

########### Regression
## Not interval censored:
survobj<-Surv(time=sero_sub_final_survival$time_obs_fail, event=sero_sub_final_survival$status2)

#make kaplan meier object
fit1<-survfit(survobj ~1,data = sero_sub_final_survival)
# fit weibull
SurvMod <- survival::survreg(survobj ~ 1,
                             dist="weibull",
                             data = sero_sub_final_survival)
summary(SurvMod)

## extract weibull params
# survreg's scale = 1/(rweibull shape)
wshape<-as.numeric(1/exp(SurvMod$icoef[2]))

# survreg's intercept = log(rweibull scale)
wscale<-exp(SurvMod$icoef[1])

## fitted 'survival'
t<-seq(0,max(sero_sub_final_survival$days_post_symptoms),0.5)
weib<-exp(-(t/wscale)^wshape)   # cumulative weibull

#str(SurvMod)
ggsurvplot(fit1,data=sero_sub_final_survival,ylab="prob still positive",xlab="days")
par(mfrow=c(1,2))
hist(rweibull(10000,shape=wshape,scale=wscale),xlab="days",main="")
plot(t,weib,type="l",xlab="days",ylab="fitted weibull curve",ylim=c(0,1))



############## With interval censoring
### don't need to specify event as we know that if time2 is missing, there was no event.
survobj<-Surv(time=sero_sub_final_survival$time1, time2=sero_sub_final_survival$time1, type = "interval2" )
fit1<-survfit(survobj ~1,data = sero_sub_final_survival)
SurvMod <- survival::survreg(survobj ~ 1,
                             dist="weibull",
                             data = sero_sub_final_survival)

summary(SurvMod)

## extract weibull params
# survreg's scale = 1/(rweibull shape)
wshape<-as.numeric(1/exp(SurvMod$icoef[2]))

# survreg's intercept = log(rweibull scale)
wscale<-exp(SurvMod$icoef[1])

## fitted 'survival'
t<-seq(0,max(sero_sub_final_survival$days_post_symptoms),0.5)
weib<-exp(-(t/wscale)^wshape)   # cumulative weibull

#str(SurvMod)
ggsurvplot(fit1,data=sero_sub_final_survival,ylab="prob still positive",xlab="days")
par(mfrow=c(1,2))
hist(rweibull(10000,shape=wshape,scale=wscale),xlab="days",main="")
plot(t,weib,type="l",xlab="days",ylab="fitted weibull curve",ylim=c(0,1))


