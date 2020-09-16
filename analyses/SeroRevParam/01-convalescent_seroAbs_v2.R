#....................................................................................................
## Purpose: Determine the distribution for the onset-of-symptoms to seroreversion
##
## Notes: Data shared from Muecksch et. al 2020
#....................................................................................................
library(tidyverse)
library(survival)
library(survminer)
source("R/my_themes.R")
source("R/extra_plotting_functions.R")
set.seed(48)

#......................
# read data
#......................
serotime <- readxl::read_excel("data/raw/shared/2020_09_11_97_for_JID.xlsx", sheet = 2) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  magrittr::set_colnames(gsub(" ", "_", colnames(.))) %>%
  magrittr::set_colnames(gsub("/", "_", colnames(.))) %>%
  magrittr::set_colnames(gsub("≥1", "gt1", colnames(.))) %>%
  dplyr::rename(sex = gender,
                donor_id = sr1407,
                days_post_symptoms = post_sx_days) %>%
  dplyr::mutate(sex = factor(sex, levels = c("F", "M")),
                donor_id = factor(donor_id))    # re-center days to months

# take to long format
serotime <- serotime %>%
  dplyr::select(-c("siemens")) %>%
  tidyr::pivot_longer(., cols = c("abbott_s_c", "diasorin_au_ml", "siemens_no.", "roche_gt1"),
                      names_to = "assay", values_to = "titres") %>%
  dplyr::mutate(assay = stringr::str_split_fixed(assay, "_", n = 2)[,1],
                assay = factor(assay, levels = c("diasorin", "siemens", "abbott", "roche"),
                               labels = c("Diasorin", "Siemens", "Abbott", "Roche")))


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
  ggplot() +
  geom_line(aes(x = post_pcr_days, y = titres, group = donor_id, color = assay)) +
  geom_hline(aes(yintercept = threshold), linetype = "dashed") +
  scale_color_manual(values = c("#3182bd", "#31a354", "#de2d26", "#756bb1")) +
  facet_wrap(~assay, scales = "free_y") +
  xyaxis_plot_theme

# look at post sx instead of pcr
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
#---- Survival Analysis #----
# Here we will subset just the abbot assay, which has the
# largest declines in sensitivity
#...........................................................
sero_final_survival <- serotime %>%
  dplyr::filter(assay == "Abbott") %>%
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(status = ifelse(titres < 1.4, 1, 0),
                status2 = ifelse(min(titres, na.rm = T) < 1.4, 1, 0),
                max_time = ifelse(days_post_symptoms == max(days_post_symptoms), 1, 0)) %>%
  dplyr::arrange(donor_id, days_post_symptoms)

# those without event. Time1 is final observation time, time2 is missing.
# not donor id 46 and 52 have missing days since symptom onset across board
sero_pos <- sero_final_survival %>%
  dplyr::filter(status2 == 0,
                !is.na(days_post_symptoms)) %>% # remove observations that weren't there
  dplyr::group_by(donor_id) %>%
  dplyr::filter(days_post_symptoms == max(days_post_symptoms)) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(time1 = days_post_symptoms,
                time2 = NA)

# those who serorevert
sero_rev <- sero_final_survival %>%
  dplyr::filter(status2 == 1) %>%
  dplyr::filter(!is.na(days_post_symptoms)) # drop observations symtpoms

# extract last time observed postiive
sero_rev_time1 <- sero_rev %>%
  dplyr::filter(status == 0) %>%
  dplyr::group_by(donor_id) %>%
  dplyr::filter(days_post_symptoms == max(days_post_symptoms)) %>%
  dplyr::mutate(time1 = days_post_symptoms)

# extract first time observed seroreverted (checked that no one becomes positive again after first being negative)
sero_rev_time2 <- sero_rev %>%
  dplyr::filter(status == 1) %>%
  dplyr::group_by(donor_id) %>%
  dplyr::filter(days_post_symptoms == min(days_post_symptoms)) %>%
  dplyr::mutate(time2=days_post_symptoms) %>%
  dplyr::select(donor_id,time2)

# combine serorev
sero_rev_comb <- dplyr::left_join(sero_rev_time1, sero_rev_time2)

# combine seropos and serorev
sero_sub_final_survival <- rbind(sero_pos, sero_rev_comb)
sero_sub_final_survival <-sero_sub_final_survival %>%
  dplyr::mutate(time_obs_fail = ifelse(status2==1, time2, time1))

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
survobj<-Surv(time=sero_sub_final_survival$time1, time2=sero_sub_final_survival$time2, type = "interval2" )
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


## Sanity check - what would the earliest possible failure time graph look like?
#Not interval censored:
survobj<-Surv(time=sero_sub_final_survival$time1, event=sero_sub_final_survival$status2)

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

