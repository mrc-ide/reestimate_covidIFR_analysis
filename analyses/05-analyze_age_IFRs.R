##................................................................................................
## Purpose: Analyze age specific IFRs
##
## Notes:
##................................................................................................
library(tidyverse)
library(COVIDCurve)
source("R/my_themes.R")
source("R/covidcurve_helper_functions.R")
source("R/delta_method.R")

# colors
study_cols <- readr::read_csv("data/plot_aesthetics/color_studyid_map.csv")
mycolors <- study_cols$cols
names(mycolors) <- study_cols$study_id
studyidnames <- study_cols %>%
  dplyr::select(c("study_id", "names")) %>%
  dplyr::filter(!is.na(names))

#......................
# read in descriptive data
#......................
dscdat <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS")

dsc_agedat <- dscdat %>%
  dplyr::filter(breakdown == "ageband") %>%
  dplyr::filter(!grepl("_nch", study_id))

#......................
# read in fitted data
#......................
regrets <- list.files("results/Modfits_noserorev/", full.names = T)
regretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(regrets), "_age", simplify = T)[,1]),
                            sero = "reg",
                            path = regrets)

serorevrets <- list.files("results/ModFits_SeroRev/", full.names = T)
serorevretmap <- tibble::tibble(study_id = toupper(stringr::str_split(basename(serorevrets), "_age", simplify = T)[,1]),
                                sero = "serorev",
                                path = serorevrets)

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
retmapSero <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(strataIFRret = purrr::map(path, get_strata_seroprevs))

#......................
# get simple seropred
#......................
simp_seroprevdat <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "seroprev_adjdat")) %>%
  tidyr::unnest(cols = "seroprev_adjdat") %>%
  dplyr::group_by(study_id, location, ageband) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # subset to latest serostudy
  group_by(study_id, location, seromidpt, ageband) %>%
  dplyr::summarise(n_totpos = sum(n_positive),
                   n_tottest = sum(n_tested)) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(sero_midday = seromidpt + lubridate::ymd("2020-01-01"),
                seroprev = n_totpos/n_tottest)
# manual liftover for studies that only reported CIs
ci_simp_seroprevdat <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "seroprev_adjdat")) %>%
  tidyr::unnest(cols = "seroprev_adjdat") %>%
  dplyr::group_by(study_id, location, ageband) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # subset to latest serostudy
  dplyr::filter(is.na(n_positive) & is.na(n_tested)) %>%
  dplyr::summarise(seroprev_ci = mean(seroprev))

# fix and take to simple column
simp_seroprevdat <- simp_seroprevdat %>%
  dplyr::left_join(., ci_simp_seroprevdat, by = c("study_id", "location", "ageband")) %>%
  dplyr::mutate(seroprev = ifelse(is.na(seroprev), seroprev_ci, seroprev)) %>%
  dplyr::select(c("study_id", "location", "ageband", "seromidpt", "sero_midday", "seroprev")) %>%
  dplyr::rename(obs_input_seropev = seroprev)

# now bring together
middays <- simp_seroprevdat %>%
  dplyr::select(c("study_id", "seromidpt")) %>%
  dplyr::filter(!duplicated(.))

# am going to need data dictionaries
datdict <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(dictkey = purrr::map(path, get_data_dict)) %>%
  tidyr::unnest(cols = "dictkey")

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
    UCI = quantile(RG_pd_seroprev, 0.975, na.rm = T),
  )  %>%
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
# Delta method needs standard error of seroprev SE(p)
# where SE(p) is the standard error of the binomial proportion for all studies except ITA
# for ITA, DNK, SWE, we use SE(p) as the logit transformed SE gleaned from the provided CIs
SE_subn <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "seroprev_adjdat")) %>%
  tidyr::unnest(cols = "seroprev_adjdat") %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::ungroup(.) %>%
  dplyr::filter(!is.na(n_positive) & !is.na(n_tested)) %>%
  dplyr::mutate(seroprev = n_positive/n_tested,
                binom_se = sqrt(seroprev * (1-seroprev))) %>%
  dplyr::select(c("study_id", "location", "ageband", "binom_se"))

SE_cis <- dsc_agedat %>%
  dplyr::select(c("study_id", "location", "seroprev_adjdat")) %>%
  tidyr::unnest(cols = "seroprev_adjdat") %>%
  dplyr::filter(study_id %in% c("ITA1", "SWE21", "DNK1")) %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::mutate(binom_se = (COVIDCurve:::logit(serouci) - COVIDCurve:::logit(serolci))/(1.96 * 2))  %>%
  dplyr::select(c("study_id", "location", "ageband", "binom_se"))

ageSE <- dplyr::bind_rows(SE_subn, SE_cis)

# get delta crude IFRs
crude_IFRs <- dsc_agedat %>%
  dplyr::select(c("study_id", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat") %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::filter(obsday == seromidpt) %>%  # sero obs day
  dplyr::ungroup(.) %>%
  dplyr::select(c("study_id", "ageband", "cumdeaths", "popn", "seroprev")) %>%
  dplyr::left_join(., ageSE, by = c("study_id", "ageband")) %>%
  dplyr::group_by_at(c("study_id", "ageband")) %>%
  dplyr::mutate(IFRcalc = cumdeaths  / (seroprev * popn + cumdeaths),
                IFRbound = purrr::map(seroprev, get_delta_CI_vals, deaths = cumdeaths, popN = popn, SE = binom_se, tol = 1e-4),
                lower_ci = purrr::map_dbl(IFRbound, "lower.ci"),
                upper_ci = purrr::map_dbl(IFRbound, "upper.ci")) %>%
  dplyr::ungroup(.) %>%
  dplyr::mutate(IFRcalc = round(IFRcalc * 100, 2),
                lower_ci = round(lower_ci * 100, 2),
                upper_ci = round(upper_ci * 100, 2))
crude_IFRs_column <- crude_IFRs %>%
  dplyr::mutate(crude_IFRs = paste0(IFRcalc, " (", lower_ci, ", ", upper_ci, ")")) %>%
  dplyr::select(c("study_id", "ageband", "crude_IFRs"))


#......................
# make final table
#......................
dplyr::left_join(simp_seroprevdat, seroprev_age_column, by = c("study_id", "ageband")) %>%
  dplyr::left_join(., crude_IFRs_column, by = c("study_id", "ageband")) %>%
  dplyr::left_join(., modIFR_age_column, by = c("study_id", "ageband")) %>%
  dplyr::select(c("location", "study_id", "sero_midday", "obs_input_seropev", "ageband", "seroprev_reg", "seroprev_serorev", "crude_IFRs", "regIFR", "serorevIFR")) %>%
  dplyr::mutate(age_low = as.numeric(stringr::str_split_fixed(ageband, "-", n=2)[,1])) %>%
  dplyr::arrange(location, study_id, age_low) %>%
  dplyr::select(-c("age_low")) %>%
  readr::write_tsv(., path = "tables/final_tables/age_specific_ifr_data.tsv")







#............................................................
#---- Log Fig Age Specific Results #----
#...........................................................

#......................
# get log transformed variables
#......................
log_retmapIFR <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(strataIFRret_log = purrr::map(path, get_log_transformed_IFR_cred_intervals, by_chain = F))

log_retmapIFR_dat <- log_retmapIFR %>%
  tidyr::unnest(cols = "strataIFRret_log") %>%
  dplyr::rename(Strata = param) %>%
  dplyr::left_join(., datdict) %>%
  dplyr::mutate(age_mid = purrr::map_dbl(ageband, function(x){
    nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
    nums[nums == 999] <- 100
    return(mean(nums))}))


#......................
# perform log-linear regression from posteriors
#......................
# NO serorev
log_regressdat_NOseroreov <- log_retmapIFR_dat %>%
  dplyr::filter(sero == "reg") %>%
  dplyr::select(c("study_id", "age_mid", "median", "precision"))

NoSerorv_model_weighted <- lm(median ~ age_mid,
                              data = log_regressdat_NOseroreov,
                              weights = precision)

# serorev
log_regressdat_seroreov <- log_retmapIFR_dat %>%
  dplyr::filter(sero == "serorev") %>%
  dplyr::select(c("study_id", "age_mid", "median", "precision"))

Serorv_model_weighted <- lm(median ~ age_mid,
                              data = log_regressdat_seroreov,
                              weights = precision)



#......................
# plot out
#......................
newdat <- data.frame(age_mid = seq(0, 95, by = 0.05))
NoSerorv_model_weighted_CIband <- predict(NoSerorv_model_weighted,
                                          newdat,
                                          interval = "confidence", alpha = 0.05)
NoSerorv_model_weighted_CIband <- cbind.data.frame(NoSerorv_model_weighted_CIband,
                                                   NoSerorv_model_weighted_CIband)

log_age_IFR_plotObj <- log_retmapIFR_dat %>%
  dplyr::filter(sero == "reg") %>%
  ggplot() +
  geom_ribbon(data = newdat, aes(x = age_mid, y = fit, ymin = lwr, ymax = upr),
              alpha = 0.35, color = "#d9d9d9", linetype = "dashed") +
  geom_pointrange(aes(x = age_mid, y = median, ymin = LCI, ymax = UCI, color =  study_id),
                  alpha = 0.75, shape = 16, size = 0.9) +
  scale_color_manual("Study ID", values = mycolors) +
  ylab("Logged Age-Specific IFR (95% CrI)") + xlab("Mid. Age") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))

#......................
# tidy modelled and plot
#......................
modIFR_age <- modIFR_age %>%
  dplyr::mutate(age_mid = purrr::map_dbl(ageband, function(x){
    nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
    nums[nums == 999] <- 100
    return(mean(nums))})) %>%
  dplyr::filter(sero == "reg")

modIFR_age_plotObj <- modIFR_age %>%
  dplyr::filter(sero == "reg") %>%
  ggplot() +
  geom_pointrange(aes(x = age_mid, y = median, ymin = LCI, ymax = UCI, color =  study_id),
                  alpha = 0.75, shape = 16, size = 0.9) +
  scale_color_manual("Study ID", values = mycolors) +
  ylab("Age-Specific IFR (95% CrI)") + xlab("Mid. Age") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))


#......................
# bring together
#......................
FigA <- modIFR_age_plotObj + theme(legend.position = "none")
FigB <- log_age_IFR_plotObj + theme(legend.position = "none")
mainFig <- cowplot::plot_grid(FigA, FigB,
                              align = "h", ncol = 2, nrow = 1,
                              labels = c("(A)", "(B)"))
legend <- cowplot::get_legend(modIFR_age_plotObj +
                                theme(legend.position = "bottom"))
mainFig <- cowplot::plot_grid(mainFig, legend,
                              nrow = 2, ncol = 1, rel_heights = c(0.9, 0.1))

jpeg("figures/final_figures/IFR_age_spec_logplot.jpg",
     width = 11, height = 8, units = "in", res = 600)
plot(mainFig)
graphics.off()




#............................................................
#---- Best IFR Est Table  #----
#...........................................................
# using log-linear regression from posteriors above

# new dat is 10 year age bands
newdat <- data.frame(age_mid = seq(5, 95, by = 10))

NoSerorv_model_weighted
Serorv_model_weighted



#......................
# bring together
#......................
NoSerorv_bestest <- predict(NoSerorv_model_weighted,
                                          newdat,
                                          interval = "confidence", alpha = 0.05) %>%
  dplyr::mutate(bestest = paste0(fit, "(", lwr, ", ", upr, ")")) %>%
  dplyr::pull(bestest)

Serorv_bestest <- predict(Serorv_model_weighted,
                            newdat,
                            interval = "confidence", alpha = 0.05) %>%
  dplyr::mutate(bestest = paste0(fit, "(", lwr, ", ", upr, ")")) %>%
  dplyr::pull(bestest)

best_est <- tibble::tibble(ageband = paste0(seq(0, 90, 10), "-", seq(10, 100, 10)),
                           NoSeroRev = NoSerorv_bestest,
                           SeroRev = Serorv_bestest)
# overall
summary(NoSerorv_model_weighted)
summary(Serorv_model_weighted)
confint(NoSerorv_model_weighted, alpha = 0.05)
confint(Serorv_model_weighted, alpha = 0.05)
