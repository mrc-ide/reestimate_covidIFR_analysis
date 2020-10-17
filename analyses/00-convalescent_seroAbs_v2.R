#....................................................................................................
## Purpose: Determine the distribution for the onset-of-symptoms to seroreversion
## using the Weibull Distribution Fit
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


neg_at_baseline_id_drop <- neg_at_baseline %>%
  dplyr::ungroup(.) %>%
  dplyr::filter(neg_at_baseline == T) %>%
  dplyr::select(c("assay", "donor_id")) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::mutate(drop = TRUE)

# do negative at baseline ever become positive? --> no
serotime %>%
  dplyr::left_join(., neg_at_baseline_id_drop, by = c("assay", "donor_id")) %>%
  dplyr::filter(!is.na(drop)) %>%
  ggplot() +
  geom_line(aes(x = days_post_symptoms, y = titres, group = donor_id, color = assay),
            color = "#3182bd") +
  geom_point(aes(x = days_post_symptoms, y = titres, group = donor_id, color = assay),
             color = "#3182bd") +
  geom_hline(yintercept = 1.4, linetype = "dashed") +
  facet_wrap(~assay, scales = "free_y") +
  ylab("Ab. Titres") + xlab("Days Post-Symptom Onset") + ggtitle("Final Included") +
  xyaxis_plot_theme

# drop those individuals negative at baseline
serotime <- serotime %>%
  dplyr::left_join(., neg_at_baseline_id_drop, by = c("assay", "donor_id")) %>%
  dplyr::filter(is.na(drop)) %>%
  dplyr::select(-c("drop"))

# quick look at finals
serotime %>%
  dplyr::left_join(., thresholds, by = "assay") %>%
  dplyr::filter(assay == "Abbott") %>%
  ggplot() +
  geom_line(aes(x = days_post_symptoms, y = titres, group = donor_id, color = assay),
            color = "#3182bd") +
  geom_point(aes(x = days_post_symptoms, y = titres, group = donor_id, color = assay),
             color = "#3182bd") +
  geom_hline(aes(yintercept = threshold), linetype = "dashed") +
  facet_wrap(~assay, scales = "free_y") +
  ylab("Ab. Titres") + xlab("Days Post-Symptom Onset") + ggtitle("Final Included") +
  xyaxis_plot_theme


#............................................................
# quick look at characteristics
#...........................................................
serotime_char <- serotime %>%
  dplyr::select(c("age", "sex", "donor_id")) %>%
  dplyr::filter(!duplicated(.))
table(serotime_char$sex)
summary(serotime_char$age)

#............................................................
#---- Survival Analysis #----
#...........................................................
# Wrangle Survival Data
# know from above no one seroconverts later than start and that no one
# who was -ve at baseline later converts

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

# get last observation date as time to event for seroNEVERts
max_time_to_event_neg <- sero_final_survival %>%
  dplyr::group_by(donor_id) %>%
  dplyr::filter(status == 0) %>%
  dplyr::filter(! donor_id %in% posdonors) %>%
  dplyr::summarise(time_to_event = max(days_post_symptoms),
                   time_to_event2 = NA) %>%
  dplyr::mutate(status = 0)



# combine serorev
sero_rev_comb <- dplyr::bind_rows(max_time_to_event_neg, min_time_to_event_pos)

# for no interval censoring account for time of failure
sero_rev_comb <- sero_rev_comb %>%
  dplyr::mutate(time_obs_fail = ifelse(status == 1, time_to_event2, time_to_event))

#............................................................
# Weibull Regression
#...........................................................
#......................
# NO interval censoring
#......................
survobj_rcens <- survival::Surv(time = sero_rev_comb$time_obs_fail,
                          event = sero_rev_comb$status)

# kaplan meier fit
KM1_mod <- survival::survfit(survobj_rcens ~ 1, data = sero_rev_comb)
summary(KM1_mod)

# fit weibull
WBmod1 <- survival::survreg(survobj_rcens ~ 1,
                           dist="weibull",
                           data = sero_rev_comb)
summary(WBmod1)


#......................
# WITH interval censoring
#......................
survobj_intcens <- survival::Surv(time = sero_rev_comb$time_to_event,
                          time2 = sero_rev_comb$time_to_event2,
                          type = "interval2" )

# fit KM
KM2_mod <- survival::survfit(survobj_intcens ~ 1,
                             data = sero_rev_comb)
summary(KM2_mod)



# fit weibull
WBmod2 <- survival::survreg(survobj_intcens ~ 1,
                           dist="weibull",
                           data = sero_rev_comb)
summary(WBmod2)


#............................................................
# Extract Weibull Params (with interval censoring)
#...........................................................
# survreg's scale = 1/(rweibull shape)
# survreg's intercept = log(rweibull scale)
weibull_params <- list(wshape = 1/exp(WBmod2$icoef[2]),
                       wscale = exp(WBmod2$icoef[1]))

# save out parameters
dir.create(path = "results/prior_inputs/", recursive = TRUE)
saveRDS(weibull_params, "results/prior_inputs/weibull_params.RDS")

#............................................................
#---- Figure of Seroreversion #----
#...........................................................
#......................
#  Kaplan Meier plot
#......................
KMplot <- survminer::ggsurvplot(fit = KM2_mod)

#......................
# Weibull plot
#......................
# fitted 'survival'
# https://stackoverflow.com/questions/9151591/how-to-plot-the-survival-curve-generated-by-survreg-package-survival-of-r
pw <- seq(from = 0, to = 1, by = 0.01)
tof_weibull <- tibble::tibble(prob = 1 - pw,
                              tof = predict(WBmod2, type="quantile", p = pw)[1,])

# KM pieces
survdat <- KMplot$data.survplot
runin <- tibble::tibble(time = c(0, min(survdat$time)), surv = c(1, 1))
survdat <- dplyr::bind_rows(runin, survdat)
censored <- KMplot$data.survplot %>%
  dplyr::filter(n.censor != 0)
events <- KMplot$data.survplot %>%
  dplyr::filter(n.event != 0)


# polotObj pieces
WeibullSurvPlotObj <- ggplot() +
  geom_line(data = survdat, aes(x = time, y = surv),
            color = "#3C3B6E", alpha = 0.9, size = 1.2) +
  geom_ribbon(data = KMplot$data.survplot, aes(x = time, ymin = lower, ymax = upper),
              fill = "#6967bf", alpha = 0.5) +
  geom_point(data = censored, aes(x = time, y = surv, size = n.censor),
             color = "#3C3B6E", alpha = 0.9, shape = 124, show.legend = F) +
  geom_point(data = events, aes(x = time, y = surv, size = n.event),
             color = "#3C3B6E", alpha = 0.9, shape = 19, show.legend = F) +
  geom_line(data = tof_weibull, aes(x = tof, y = prob),
            size = 1.25, color = "#B22234", alpha = 0.8) +
  ylab("Prob. of Seropositivity Persistence") +
  xlab("Time (Days)") +
  xyaxis_plot_theme


dir.create("figures/final_figures/", recursive = TRUE)
jpeg("figures/final_figures/weibull_survplot.jpg",
     width = 11, height = 8, units = "in", res = 500)
plot(WeibullSurvPlotObj)
graphics.off()


# save out
saveRDS(WeibullSurvPlotObj, "figures/final_figures/weibull_survplot.RDS")

