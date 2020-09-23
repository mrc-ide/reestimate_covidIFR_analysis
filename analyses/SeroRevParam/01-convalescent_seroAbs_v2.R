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
# Here we will subset just the abbot assay, which has the
# largest declines in sensitivity
# also will drop missing post-sx observation times
#...........................................................
# three missing post day sxs
serotime <- serotime %>%
  dplyr::filter(!is.na(days_post_symptoms)) %>%
  dplyr::filter(assay == "Abbott")



#............................................................
# Exclude individuals negative at baseline
#...........................................................
neg_at_baseline <-  serotime %>%
  dplyr::filter(!is.na(days_post_symptoms)) %>% # remove missing post symp days
  dplyr::group_by(donor_id, assay) %>%
  dplyr::filter(days_post_symptoms == min(days_post_symptoms)) %>%
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
serotime <- serotime %>%
  dplyr::left_join(., neg_at_baseline, by = c("assay", "donor_id")) %>%
  dplyr::filter(is.na(drop)) %>%
  dplyr::select(-c("drop"))

# quick look at finals
serotime %>%
  dplyr::left_join(., thresholds, by = "assay") %>%
  dplyr::filter(assay == "Abbott") %>%
  ggplot() +
  geom_line(aes(x = days_post_symptoms, y = titres, group = donor_id, color = assay),
            color = "#3182bd") +
  geom_hline(aes(yintercept = threshold), linetype = "dashed") +
  facet_wrap(~assay, scales = "free_y") +
  ylab("Ab. Titres") + xlab("Days Post-Symptom Onset") + ggtitle("Final Included") +
  xyaxis_plot_theme


#............................................................
#---- Survival Analysis #----
#...........................................................
#............................................................
# Wrangle Survival Data
#...........................................................
sero_final_survival <- serotime %>%
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(status = ifelse(titres < 1.4, 1, 0)) %>%
  dplyr::arrange(donor_id, days_post_symptoms) %>%
  dplyr::ungroup(.)

# get min time to event for positives
posdonors <- sero_final_survival %>%
  dplyr::filter(status == 1) %>%
  dplyr::select(c("donor_id")) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::pull("donor_id")

# time 1 for interval
min_righttime_to_event_pos <- sero_final_survival %>%
  dplyr::filter(donor_id %in% posdonors) %>%
  dplyr::filter(status == 1) %>%
  dplyr::group_by(donor_id) %>%
  dplyr::summarise(time_to_event2 = min(days_post_symptoms)) %>%
  dplyr::mutate(status = 1)

# time w for interval
min_lefttime_to_event_pos <- sero_final_survival %>%
  dplyr::filter(donor_id %in% posdonors) %>%
  dplyr::group_by(donor_id) %>%
  dplyr::filter(status == 0) %>%
  dplyr::summarise(time_to_event = max(days_post_symptoms)) %>%
  dplyr::select(c("donor_id", "time_to_event"))


min_time_to_event_pos <- dplyr::left_join(min_righttime_to_event_pos,
                                          min_lefttime_to_event_pos, by = "donor_id")

# get last observation date as time to event for never seroreverteds
max_time_to_event_neg <- sero_final_survival %>%
  dplyr::group_by(donor_id) %>%
  dplyr::filter(status == 0) %>%
  dplyr::filter(! donor_id %in% posdonors) %>%
  dplyr::summarise(time_to_event = max(days_post_symptoms),
                   time_to_event2 = NA) %>%
  dplyr::mutate(status = 0)

#......................
# exclude individuals who serorevert but then become seropositive later
#......................
excl_inds <- sero_final_survival %>%
  dplyr::filter(titres < 1.14) %>%
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(lag_time = days_post_symptoms - dplyr::lag(days_post_symptoms))
# looks good to go


# combine serorev
sero_rev_comb <- dplyr::bind_rows(max_time_to_event_neg, min_time_to_event_pos)

# for no interval censoring account for time of faiure
sero_rev_comb <- sero_rev_comb %>%
  dplyr::mutate(time_obs_fail = ifelse(status == 1, time_to_event2, time_to_event))

#............................................................
# Weibull Regression
#...........................................................
#......................
# NO interval censoring
#......................
survobj_km <- survival::Surv(time = sero_rev_comb$time_obs_fail,
                          event = sero_rev_comb$status)

#make kaplan meier object
KM1_mod <- survival::survfit(survobj_km ~ 1, data = sero_rev_comb)

# fit weibull
WBmod <- survival::survreg(survobj_km ~ 1,
                           dist="weibull",
                           data = sero_rev_comb)
summary(WBmod)


#......................
# WITH interval censoring
#......................
survobj <- survival::Surv(time = sero_rev_comb$time_to_event,
                          time2 =sero_rev_comb$time_to_event2,
                          type = "interval2" )

# fit weibull
WBmod <- survival::survreg(survobj ~ 1,
                           dist="weibull",
                           data = sero_rev_comb)
summary(WBmod)


#............................................................
# Extract Weibull Params (with interval censoring)
#...........................................................
# survreg's scale = 1/(rweibull shape)
# survreg's intercept = log(rweibull scale)
weibull_params <- list(wshape = 1/exp(WBmod$icoef[2]),
                       wscale = exp(WBmod$icoef[1]))

#......................
# save out results
#......................
dir.create(path = "results/sero_reversion/", recursive = TRUE)
saveRDS(sero_rev_comb, "results/sero_reversion/sero_rev_dat.RDS")
saveRDS(weibull_params, "results/sero_reversion/weibull_params.RDS")
# NB -- KM represents the earliest possible failure time graph
saveRDS(KM1_mod, "results/sero_reversion/KaplanMeierFit.RDS")
saveRDS(survobj_km, "results/sero_reversion/survobj_km.RDS")
saveRDS(WBmod, "results/sero_reversion/WeibullFit.RDS")
saveRDS(serotime, file = "results/sero_reversion/sero_reversion_incld_data.RDS")

