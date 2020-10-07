## ....................................................................................................
## Purpose:
## Reads in seroconversion/reversion data and fits the mean time to reversion
## (mu) under a simple model of constant hazard lambda of seroconverting and
## constant hazard mu of seroreverting.
##
## Notes: Data shared from Muecksch et. al 2020
## ....................................................................................................
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
#...........................................................
serotime <- serotime %>%
  dplyr::filter(!is.na(post_pcr_days)) %>%
  dplyr::filter(assay == "Abbott")



#............................................................
# Exclude individuals negative at baseline
#...........................................................
neg_at_baseline <-  serotime %>%
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
  geom_line(aes(x = post_pcr_days, y = titres, group = donor_id, color = assay),
            color = "#3182bd") +
  geom_point(aes(x = post_pcr_days, y = titres, group = donor_id, color = assay),
            color = "#3182bd") +
  geom_hline(aes(yintercept = threshold), linetype = "dashed") +
  facet_wrap(~assay, scales = "free_y") +
  ylab("Ab. Titres") + xlab("Days Post-PCR") + ggtitle("Final Included") +
  xyaxis_plot_theme



#......................
# does anyone serorevert but then become seropositive later --> no
#......................
serorevert_ids <-  serotime %>%
  dplyr::group_by(donor_id, assay) %>%
  dplyr::summarise(mintitre = min(titres)) %>%
  dplyr::filter(mintitre <= 1.4) %>%
  dplyr::pull(donor_id)
serorevert_ids <- as.character( serorevert_ids[!duplicated(serorevert_ids)] )
serotime %>%
  dplyr::filter(donor_id %in% serorevert_ids) %>%
  ggplot() +
  geom_line(aes(x = days_post_symptoms, y = titres, group = donor_id, color = assay),
            color = "#3182bd") +
  geom_point(aes(x = days_post_symptoms, y = titres, group = donor_id, color = assay),
             color = "#3182bd") +
  geom_hline(yintercept = 1.4, linetype = "dashed") +
  facet_wrap(~assay, scales = "free_y") +
  ylab("Ab. Titres") + xlab("Days Post-Symptom Onset") +
  xyaxis_plot_theme

#............................................................
# quick look at final group characteristics
#...........................................................
serotime_char <- serotime %>%
  dplyr::select(c("age", "sex", "donor_id")) %>%
  dplyr::filter(!duplicated(.))
table(serotime_char$sex)
summary(serotime_char$age)

#............................................................
# redistribute time from seroconversion to serorev
#...........................................................
# for each individual, find time of seroconversion
# know from figure, no one seroconverts then seroreverts then seroconverts
serocon <- serotime %>%
  dplyr::group_by(donor_id) %>%
  dplyr::filter(titres >= 1.4) %>%
  dplyr::summarise(day_of_serocon = min(post_pcr_days))

# add in
serotime <- serotime %>%
  dplyr::left_join(., serocon, by = "donor_id") %>%
  dplyr::mutate(days_since_serocon = post_pcr_days - day_of_serocon)

# plot out
serotime %>%
  ggplot() +
  geom_point(aes(x = days_since_serocon, y = day_of_serocon))

# tidy up
sero_final_survival <- serotime %>%
  dplyr::group_by(donor_id) %>%
  dplyr::mutate(status = ifelse(titres < 1.4, 1, 0)) %>%
  dplyr::arrange(donor_id, days_since_serocon) %>%
  dplyr::ungroup(.)

# get min time to event for positives
posdonors <- sero_final_survival %>%
  dplyr::filter(status == 1) %>%
  dplyr::select(c("donor_id")) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::pull("donor_id")

# time 2 for interval
min_righttime_to_event_pos <- sero_final_survival %>%
  dplyr::filter(donor_id %in% posdonors) %>%
  dplyr::filter(status == 1) %>%
  dplyr::group_by(donor_id) %>%
  dplyr::summarise(time_to_event2 = min(days_since_serocon)) %>%
  dplyr::mutate(status = 1)

# time 1 for interval
min_lefttime_to_event_pos <- sero_final_survival %>%
  dplyr::filter(donor_id %in% posdonors) %>%
  dplyr::group_by(donor_id) %>%
  dplyr::filter(status == 0) %>%
  dplyr::summarise(time_to_event = max(days_since_serocon)) %>%
  dplyr::select(c("donor_id", "time_to_event"))


min_time_to_event_pos <- dplyr::left_join(min_righttime_to_event_pos,
                                          min_lefttime_to_event_pos, by = "donor_id")

# get last observation date as time to event for never seroreverteds
max_time_to_event_neg <- sero_final_survival %>%
  dplyr::group_by(donor_id) %>%
  dplyr::filter(! donor_id %in% posdonors) %>%
  dplyr::summarise(time_to_event = NA,
                   time_to_event2 = max(days_since_serocon)) %>%
  dplyr::mutate(status = 0)

# combine serorev
sero_rev_comb <- dplyr::bind_rows(max_time_to_event_neg, min_time_to_event_pos)



#............................................................
# Fit Seroreversion
#...........................................................
# define fixed parameters -- this is onset of infxn to seroconversion
lambda <- 13.3 + 5

# pdf of seroconverting at time t
ft <- function(t, mu) {
  1/(mu - lambda)*exp(-t/mu) - 1/(mu - lambda)*exp(-t/lambda)
}

# cdf of seroconverting at time t
Ft <- function(t, mu) {
  lambda/(mu - lambda)*(exp(-t/lambda) - 1) - mu/(mu - lambda)*(exp(-t/mu) - 1)
}

# probability of seroconverting between times t1 and t2, returned in log space
loglike <- function(t1, t2, mu) {
  ret <- Ft(t2, mu) - Ft(t1, mu)
  return(log(ret))
}

# calculate loglikelihood over range of mu
mu_vec <- seq(20, 400)
ll <- rep(0, length(mu_vec))
for (j in seq_along(mu_vec)) {
  for (i in seq_len(nrow(sero_rev_comb))) {
    t1 <- sero_rev_comb$time_to_event[i]
    t2 <- sero_rev_comb$time_to_event2[i]
    if (sero_rev_comb$status[i] == 1) {
      ll[j] <- ll[j] + loglike(t1, t2, mu_vec[j])
    } else {
      ll[j] <- ll[j] + loglike(t2, Inf, mu_vec[j])
    }
  }
}

# plot likelihood curve
plot(mu_vec, exp(ll), type = 'l')

# get maximum likelihood mu
mu_best <- mu_vec[which.max(ll)]

# save out fitted rate of seroreversion parameter
saveRDS(mu_best, "results/prior_inputs/serorev_param.RDS")


#............................................................
# Plot Out
#...........................................................
# get survival curve from complementary cdf
t <- seq(0, 300)
FSurv <- 1 - Ft(t, mu_best)
# plot Survivla curve in base
plot(t, FSurv, type = "l")
tof_fit <- data.frame(tof = t,
                      prob = 1 - Ft(t, mu_best))

#......................
# make kaplan meier object
#......................
survobj_km <- survival::Surv(time = sero_rev_comb$time_to_event2,
                             event = sero_rev_comb$status)
KM1_mod <- survival::survfit(survobj_km ~ 1, data = sero_rev_comb)

# KM plot
KMplot <- survminer::ggsurvplot(fit = KM1_mod)
# KM pieces
censored <- KMplot$data.survplot %>%
  dplyr::filter(n.censor != 0)
events <- KMplot$data.survplot %>%
  dplyr::filter(n.event != 0)
# polotObj pieces
SurvPlotObj <- ggplot() +
  geom_line(data = KMplot$data.survplot, aes(x = time, y = surv),
            color = "#3C3B6E", alpha = 0.9, size = 1.2) +
  geom_ribbon(data = KMplot$data.survplot, aes(x = time, ymin = lower, ymax = upper),
              fill = "#6967bf", alpha = 0.5) +
  geom_point(data = censored, aes(x = time, y = surv, size = n.censor),
             color = "#3C3B6E", alpha = 0.9, shape = 124, show.legend = F) +
  geom_point(data = events, aes(x = time, y = surv, size = n.event),
             color = "#3C3B6E", alpha = 0.9, shape = 19, show.legend = F) +
  geom_line(data = tof_fit, aes(x = tof, y = prob),
            size = 1.25, color = "#B22234", alpha = 0.8) +
  ylab("Prob. of Seroreversion after Seroconversion") +
  xlab("Time (Days)") +
  xyaxis_plot_theme


jpeg("figures/final_figures/exp_survplot.jpg",
     width = 11, height = 8, units = "in", res = 500)
plot(SurvPlotObj)
graphics.off()

# save out
saveRDS(SurvPlotObj, "figures/final_figures/exp_survplot.RDS")
