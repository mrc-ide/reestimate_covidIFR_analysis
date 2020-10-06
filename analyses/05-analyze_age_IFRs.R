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
#---- Read in Fitted data #----
#...........................................................
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
  dplyr::mutate(sero_midday = seromidpt + lubridate::ymd("2020-01-01") - 1,
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
source("R/monte_carlo_cis.R")
# read in observed data
dscdat <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS")
dsc_agedat <- dscdat %>%
  dplyr::filter(breakdown == "ageband") %>%
  dplyr::filter(!grepl("_nch", study_id))

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
  dplyr::mutate(age_low = as.numeric(stringr::str_split_fixed(ageband, "-", n=2)[,1])) %>%
  dplyr::mutate(obs_input_seropev = round(obs_input_seropev, 2)) %>%
  dplyr::arrange(location, study_id, age_low) %>%
  dplyr::select(-c("age_low")) %>%
  dplyr::select("location", "ageband", "sero_midday", "obs_input_seropev", "seroprev_reg", "seroprev_serorev", "crude_IFRs", "regIFR", "serorevIFR") %>%  # fix order
  readr::write_tsv(., path = "tables/final_tables/age_specific_ifr_data.tsv")


#............................................................
#---- Log Fig Age Specific Results #----
#...........................................................

#......................
# get log transformed variables
#......................
log_retmapIFR <- dplyr::bind_rows(regretmap, serorevretmap) %>%
  dplyr::mutate(strataIFRret_log = purrr::map(path, get_log10_transformed_IFR_cred_intervals, by_chain = F))

log_retmapIFR_dat <- log_retmapIFR %>%
  tidyr::unnest(cols = "strataIFRret_log") %>%
  dplyr::rename(Strata = param) %>%
  dplyr::left_join(., datdict) %>%
  dplyr::mutate(age_mid = purrr::map_dbl(ageband, function(x){
    nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
    nums[nums == 999] <- 100
    return(mean(nums))}))

#......................
# plot out
#......................

log_age_IFR_plotObj <- log_retmapIFR_dat %>%
  dplyr::filter(sero == "reg") %>%
  ggplot() +
  geom_pointrange(aes(x = age_mid, y = median, ymin = LCI, ymax = UCI, color =  location),
                  alpha = 0.75, shape = 16, size = 0.9) +
  scale_color_manual("Location", values = mycolors) +
  ylab("Log-10 Age-Specific IFR (95% CrI)") + xlab("Mid. Age") +
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
                                guides(color = guide_legend(nrow = 2)) +
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
#......................
# get "precision" based on overall IFR
#......................
studyprecision <- retmapIFR %>%
  dplyr::mutate(overallIFRret = purrr::map(path, get_overall_IFRs)) %>%
  tidyr::unnest(cols = "overallIFRret") %>%
  dplyr::mutate(wi = (UCI-LCI)) %>%
  dplyr::select(c("study_id", "wi"))

#......................
# unpack data from NLLS
#......................
ifrdat <- retmapIFR %>%
  dplyr::select(c("study_id", "sero", "strataIFRret")) %>%
  tidyr::unnest(cols = "strataIFRret") %>%
  dplyr::rename(ageband = strata) %>%
  dplyr::mutate(age_mid = purrr::map_dbl(ageband, function(x){
    nums <- as.numeric(stringr::str_split_fixed(x, "-", n = 2))
    nums[nums == 999] <- 100
    return(mean(nums))})) %>%
  dplyr::left_join(., studyprecision, by = "study_id") %>%
  dplyr::rename(mas = median) %>%
  dplyr::select(c("study_id", "sero", "ageband", "age_mid", "mas", "wi"))

#......................
# optim approach
#......................
find_pars <- function(dat) {
  if( !all(c("wi", "mas", "age_mid") %in% colnames(dat)) ) {
    stop()
  }
  # optimize exponential distribution by weighted residual sum of squares
  wi_rss <- function(par, dat) {
    sum( dat$wi * ((par[1] * exp(dat$age_mid * par[1])) - dat$mas)^2 )
  }
  # run optim
  out <- optim(par = c(1e-18), fn = wi_rss, dat = dat,
               method = "L-BFGS-B", lower = .Machine$double.xmin)
  names(out$par) <- c("lambda")
  return(out)
}

retoptim <- find_pars(ifrdat)
retoptim$par
retoptim$convergence

#......................
# get new observations
#......................
run_exp_growth <- function(par1, age){ par1 * exp(age * par1) }
newagedat <- seq(5, 95, 10)
best_est <- tibble::tibble(agemid = newagedat,
                           IFR = sapply(newagedat, run_exp_growth, par1 = retoptim$par[[1]])
                           ) %>%
  dplyr::mutate(IFR = IFR * 100)



modIFR_age_plotObj +
  geom_line(data = best_est, aes(x = agemid, y = IFR), color = "#FF0018",
            size = 1.2,
            linetype = "dashed")



# out
best_est %>%
  readr::write_tsv(., path = "tables/final_tables/overall_best_est_for_IFRs.tsv")
