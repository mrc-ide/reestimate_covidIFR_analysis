#................................................................................................
## Purpose: Plot descriptive statistics
##
## Notes:
#................................................................................................
#......................
# setup
#......................
library(tidyverse)
source("R/crude_plot_summ.R")
source("R/my_themes.R")
source("R/extra_plotting_functions.R")
dir.create("figures/descriptive_figures/", recursive = TRUE)

write2file <- F

#............................................................
#---- Read in and Wrangle Data #----
#...........................................................
# colors
study_cols <- readr::read_csv("data/plot_aesthetics/color_studyid_map.csv")
mycolors <- study_cols$cols
names(mycolors) <- study_cols$study_id

# care homes
deaths_ch <- readr::read_csv("data/raw/care_home_deaths.csv")

# data map
datmap <- readxl::read_excel("data/derived/derived_data_map.xlsx")
datmap <- datmap %>%
  dplyr::mutate(data = purrr::map(relpath, readRDS))

#......................
# wrangle & extract sero data
#......................
serohlp <- datmap %>%
  dplyr::mutate(
    seroprevdat = purrr::map(data, "seroprevMCMC"),
    sens = purrr::map(data, "sero_sens"),
    sens = purrr::map_dbl(sens, function(x){as.numeric(x$sensitivity)}),
    spec = purrr::map(data, "sero_spec"),
    spec = purrr::map_dbl(spec, function(x){as.numeric(x$specificity)})) %>%
  dplyr::select(c("seroprevdat", "sens", "spec"))

datmap <- datmap %>%
  dplyr::mutate(seroprev_adjdat = purrr::pmap(serohlp, adjust_seroprev))

#......................
# wrangle & extract death data
#......................
deathhlp <- datmap %>%
  dplyr::mutate(
    deathdat_long = purrr::map(data, "deaths_group"),
    popdat = purrr::map(data, "prop_pop"),
    groupingvar = breakdown,
    Nstandardization = 1e6) %>%
  dplyr::select(c("deathdat_long", "popdat", "groupingvar", "Nstandardization"))

datmap <- datmap %>%
  dplyr::mutate(std_deaths = purrr::pmap(deathhlp, standardize_deathdat))

#......................
# combine
#......................
datmap <- datmap %>%
  dplyr::mutate(plotdat = purrr::map2(.x = std_deaths, .y = seroprev_adjdat, dplyr::left_join)) # let dplyr find strata

# save out
dir.create("results/descriptive_results/", recursive = TRUE)
saveRDS(datmap, file = "results/descriptive_results/descriptive_results_datamap.RDS")


#............................................................
#---- Age Bands Plots/Descriptions #----
#...........................................................
#...........................................................
# GBR4 Jersey, GBR2 Scotland, SF_CA - no age data or too early in epidemic to have age data.
# LA_CA1, GBR2 TODO
#...........................................................
ageplotdat <- datmap %>%
  dplyr::filter(breakdown == "ageband") %>%
  dplyr::select(c("study_id","care_home_deaths", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat")
ageplotdat <- dplyr::full_join(ageplotdat, study_cols, by="study_id")

#filter to only plot the latest serology when there are multiple rouunds
maxDays <- ageplotdat %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(max_day=max(obsdaymax))
ageplotdat <- dplyr::full_join(ageplotdat, maxDays, by="study_id")
ageplotdat <- dplyr::filter(ageplotdat, obsdaymax == max_day)


#......................
# age raw seroprevalence
#......................
age_seroplot <- ageplotdat %>%
  dplyr::filter(care_home_deaths=="yes" & study_id!="CHE2") %>%
  dplyr::select(c("study_id", "age_mid", "seroprev")) %>%
  dplyr::mutate(seroprev = seroprev * 100) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = age_mid, y = seroprev, fill = study_id), shape = 21, size = 2.5, stroke = 0.2) +
  geom_line(aes(x = age_mid, y = seroprev, group=study_id,color=study_id), size = 0.3) +
  scale_fill_manual(values = mycolors, name = "study_id") +
  scale_color_manual(values = mycolors, name = "study_id") +
  xlab("Age (yrs).") + ylab("Raw Seroprevalence (%)") +
  xyaxis_plot_theme
if(write2file) ggsave(filename = "results/descriptive_figures/age_raw_seroplot.tiff",
                      plot = age_seroplot, width = 7, height = 5)

#......................
# age adj seroprevalence
#......................
age_seroplot <- ageplotdat %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  dplyr::select(c("study_id", "age_mid", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_line(aes(x = age_mid, y = seroprevadj, color = study_id), alpha = 0.8, size = 1.2) +
  geom_point(aes(x = age_mid, y = seroprevadj, color = study_id)) +
  scale_color_manual("Study ID", values = mycolors) +
  xlab("Age (yrs).") + ylab("Adj. Seroprevalence (%)") +
  xyaxis_plot_theme
if(write2file) ggsave(filename = "results/descriptive_figures/age_adj_seroplot.tiff",
                      plot = age_seroplot, width = 7, height = 5)

#......................
# crude IFR
#......................
# raw serology
age_IFRraw_plot0 <- ageplotdat %>%
  dplyr::filter(seromidpt == obsday) %>%
  dplyr::select(c("study_id","n_positive","n_tested","ageband", "age_mid", "cumdeaths", "popn", "seroprev", "seroprevadj","care_home_deaths")) %>%
  dplyr::mutate(infxns = popn * seroprev,
                crudeIFR =  cumdeaths/(infxns+cumdeaths),
                crudeIFR = ifelse(crudeIFR > 1, 1, crudeIFR),
                seroprev = seroprev * 100,
                sero_adj_infxns = popn *seroprevadj,
                sero_adjIFR = cumdeaths/(sero_adj_infxns+cumdeaths))

# write this out for later use
readr::write_csv(age_IFRraw_plot0, path = "data/derived/age_summ_IFR.csv")

age_IFRraw_plot <- age_IFRraw_plot0 %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  ggplot() +
  geom_line(aes(x = age_mid, y = crudeIFR, color = study_id), alpha = 0.8, size = 1.2) +
  geom_point(aes(x = age_mid, y = crudeIFR, fill = seroprev), color = "#000000", size = 2.5, shape = 21, alpha = 0.8) +
  scale_color_manual("Study ID", values = mycolors) +
  scale_fill_gradientn("Raw Seroprev.",
                       colors = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))) +
  xlab("Age (yrs).") + ylab("Crude Infection Fatality Rate") +
  xyaxis_plot_theme
if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/age_IFRraw_plot.jpg",
                           plot = age_IFRraw_plot)


### probably a neater way to do this.

age_IFRraw_plot2 <- age_IFRraw_plot0 %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = age_mid, y = crudeIFR, fill = study_id), shape = 21, size = 2.5, stroke = 0.2) +
  geom_line(aes(x = age_mid, y = crudeIFR, group=study_id,color=study_id), size = 0.3) +
  scale_fill_manual(values = mycolors, name = "study_id") +
  scale_color_manual(values = mycolors, name = "study_id") +
  xlab("Age (years)") + ylab("Crude infection fatality rate") +
  xyaxis_plot_theme
if(write2file) ggsave(filename = "results/descriptive_figures/age_IFRraw_plot2.pdf", plot = age_IFRraw_plot2, width = 7, height = 5)

age_IFRraw_plot_log <- age_IFRraw_plot0 %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  filter(crudeIFR>0) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = age_mid, y = crudeIFR, fill = study_id), shape = 21, size = 2.5, stroke = 0.2) +
  geom_line(aes(x = age_mid, y = crudeIFR, group=study_id,color=study_id), size = 0.3) +
  scale_fill_manual(values = col_vec, name = "study_id") +
  scale_color_manual(values = mycolors, name = "study_id") +
  xlab("Age (years)") + ylab("Crude infection fatality rate") +
  xyaxis_plot_theme +
  scale_y_log10()
#  coord_cartesian(ylim=c(0.00000000001,1))
if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/age_IFRraw_plot_log.jpg",
                           plot = age_IFRraw_plot_log,width_wide = 8,height_wide = 5.5)

# seroadj
age_IFRadj_plot <- age_IFRraw_plot0 %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = age_mid, y = sero_adjIFR, fill = study_id), shape = 21, size = 2.5, stroke = 0.2) +
  geom_line(aes(x = age_mid, y = sero_adjIFR, group=study_id,color=study_id), size = 0.3) +
  scale_fill_manual(values = col_vec, name = "study_id") +
  scale_color_manual(values = col_vec, name = "study_id") +
  xlab("Age (years)") + ylab("Adjusted infection fatality rate") +
  xyaxis_plot_theme
if(write2file) ggsave(filename = "results/descriptive_figures/age_IFRadj_plot.pdf", plot = age_IFRadj_plot, width = 7, height = 5)


#......................
# compare with and without care home deaths
#......................
study_ids_ch<-c(deaths_ch$study_id,paste0(deaths_ch$study_id,"_nch"))
study_cols_ch<-filter(study_cols,study_id %in% study_ids_ch)
col_vec<-study_cols_ch$study_cols
names(col_vec) <- study_cols_ch$country
age_IFRraw_plot_ch<-dplyr::filter(age_IFRraw_plot0,study_id %in% study_ids_ch & care_home_deaths=="yes")
age_IFRraw_plot_noch<-dplyr::filter(age_IFRraw_plot0,study_id %in% study_ids_ch & care_home_deaths=="no")

age_IFRraw_plot_ch <- ggplot() + theme_bw() +
  geom_point(aes(x = age_IFRraw_plot_noch$age_mid, y = age_IFRraw_plot_noch$crudeIFR, fill = age_IFRraw_plot_noch$study_id), shape = 21, size = 2.5, stroke = 0.2) +
  geom_line(aes(x = age_IFRraw_plot_noch$age_mid, y = age_IFRraw_plot_noch$crudeIFR, group=age_IFRraw_plot_noch$study_id,color=age_IFRraw_plot_noch$study_id), size = 0.3) +
  #scale_fill_manual(values = col_vec, name = "country") +
  #scale_color_manual(values = col_vec, name = "country") +
  xlab("Age (years)") + ylab("Crude infection fatality rate") +
  xyaxis_plot_theme
if(write2file) ggsave(filename = "results/descriptive_figures/age_IFRraw_plot_ch.tiff", plot = age_IFRraw_plot_ch, width = 7, height = 5)


#......................
# standardized deaths by age
#......................
age_std_cum_deaths_plot <- ageplotdat %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  dplyr::filter(seromidpt == obsday) %>%
  dplyr::select(c("study_id", "age_mid", "std_cum_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_line(aes(x = age_mid, y = std_cum_deaths, color = study_id), alpha = 0.8, size = 1.2) +
  geom_point(aes(x = age_mid, y = std_cum_deaths, fill = seroprevadj), color = "#000000", size = 2.5, shape = 21, alpha = 0.8) +
  scale_color_manual("Study ID", values = study_id) +
  scale_fill_gradientn("Adj. Seroprevalence (%)",
                       colors = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))) +
  xlab("Age (yrs).") + ylab("Cum. Deaths per Million") +
  labs(caption = "Cumulative Deaths per Million at midpoint of Seroprevalence Study") +
  xyaxis_plot_theme
if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/age_std_cum_deaths_plot.jpg",
                           plot = age_std_cum_deaths_plot)


#......................
# cumulative proportion deaths by age
#......................
age_IFRraw_plot0 <- age_IFRraw_plot0 %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  dplyr::mutate(d_per_mill=cumdeaths/popn)

tot_deaths <- age_IFRraw_plot0 %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(tot_deaths=sum(cumdeaths),
                   tot_deaths_std=sum(d_per_mill)) %>%
  dplyr::select(study_id, tot_deaths,tot_deaths_std) %>%
  ungroup()

age_prop_deaths_plotdat <- full_join(age_IFRraw_plot0,tot_deaths,by="study_id") %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  dplyr::mutate(prop_deaths = cumdeaths/tot_deaths,
                prop_deaths_std = d_per_mill/tot_deaths_std) %>%
  dplyr::arrange(study_id,age_mid)

cumu_deaths <- age_prop_deaths_plotdat %>%
  dplyr::group_by(study_id,age_mid) %>%
  dplyr::summarise(cum_prop_deaths=cumsum(prop_deaths),
                   cum_prop_deaths_std=cumsum(prop_deaths_std))

age_prop_deaths_plotdat <-age_prop_deaths_plotdat %>%
  dplyr::mutate(cum_prop_deaths=cumu_deaths$cum_prop_deaths,
                cum_prop_deaths_std=cumu_deaths$cum_prop_deaths_std,
                age_low = as.numeric(stringr::str_extract(ageband, "[0-9]+?(?=-)")),
                age_high = as.numeric(gsub( "(.*)-(.*)", "\\2",  ageband)))

prop_deaths_70<-age_prop_deaths_plotdat %>%
  dplyr::mutate(ageband2=ifelse(age_low<69,"0-69","70+")) %>%
  dplyr::group_by(study_id,ageband2) %>%
  dplyr::summarise(prop_deaths=sum(prop_deaths),
                   prop_deaths_std=sum(prop_deaths_std)) %>%
  dplyr::filter(ageband2=="70+")

########## Raw deaths by age cumulative (just showing us the population structure more than anything?)
age_prop_deaths_plot<-ggplot(age_prop_deaths_plotdat, aes(x = age_mid, y = cum_prop_deaths, group=study_id)) +
  #geom_point(aes(fill = study_id), color = "#000000", size = 2.5, shape = 21, alpha = 0.8) +
  geom_line(aes(color = study_id), alpha = 0.8, size = 1.2) +
  scale_color_manual("Study ID", values = c("ITA1"=discrete_colors[1],
                                            "ESP1-2"=discrete_colors[2],
                                            "GBR3"=discrete_colors[3],
                                            "NLD1"=discrete_colors[4],
                                            "CHN1"=discrete_colors[5],
                                            "NYC_NY_1"=discrete_colors[6],
                                            "BRA1"=discrete_colors[7],
                                            "CHE1"=discrete_colors[8],
                                            "CHE2"=discrete_colors[9],
                                            "DNK1"=discrete_colors[10],
                                            "LUX1"=discrete_colors[11])) +
  xlab("Age (yrs)") + ylab("Cumulative proportion of deaths") +
  xyaxis_plot_theme
if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/age_prop_cum_deaths_plot.jpg",
                           plot = age_prop_deaths_plot)

########## Deaths per capita by age cumulative
col_vec<-study_cols$study_cols
names(col_vec) <- study_cols$study_id
age_prop_deaths_plot_std <- age_prop_deaths_plotdat %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  ggplot() + theme_bw() +
  geom_line(aes(x = age_mid, y = cum_prop_deaths_std, group=study_id,color=study_id), size = 0.3) +
  scale_color_manual(values = col_vec, name = "study_id") +
  xlab("Age (yrs)") + ylab("Cumulative proportion of deaths, age-standardised") +
  xyaxis_plot_theme
if(write2file) ggsave(filename = "results/descriptive_figures/age_prop_deaths_plot_std.tiff", plot = age_prop_deaths_plot_std, width = 7, height = 5)

#......................
# daily standardized deaths by age
#......................
age_std_daily_deaths_plot <- ageplotdat %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  dplyr::select(c("study_id", "obsday", "ageband", "age_mid", "std_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(ageband = forcats::fct_reorder(ageband, age_mid),
                seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_line(aes(x = obsday, y = std_deaths, color = ageband), alpha = 0.8, size = 1.2) +
  facet_wrap(.~study_id, scales = "free_y") +
  xlab("Obs. Day") + ylab("Daily Deaths per Million") +
  xyaxis_plot_theme
if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/age_std_daily_deaths_plot.jpg",
                           plot = age_std_daily_deaths_plot)

#............................................................
#----  Regional Plots/Descriptions #----
#...........................................................
rgnplotdat <- datmap %>%
  dplyr::filter(breakdown == "region" & care_home_deaths=="yes") %>%
  dplyr::select(c("study_id", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat")

# filter to only plot the latest serology when there are multiple rounds
maxDays <- rgnplotdat %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(max_day=max(obsdaymax))
rgnplotdat <- dplyr::full_join(rgnplotdat,maxDays,by="study_id")
rgnplotdat <- dplyr::filter(rgnplotdat,obsdaymax==max_day)
rgnplotdat<-full_join(rgnplotdat,study_cols,by="study_id")

#......................
# rgn adj seroprevalence
#......................
# col_vec<-study_cols$study_cols
# names(col_vec) <- study_cols$study_id
rgn_seroplot <- rgnplotdat %>%
  dplyr::select(c("study_id", "region", "seroprev")) %>%
  dplyr::mutate(seroprev = seroprev * 100) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = region, y = seroprev, color = study_id), size = 2.5) +
  scale_color_manual(values = col_vec, name = "study_id") +
  facet_wrap(.~study_id, scales = "free_x") +
  xlab("Region") + ylab("Raw Seroprevalence (%)") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/rgn_raw_seroplot.jpg",
                           plot = rgn_seroplot)
########## Deaths per capita by age cumulative
age_prop_deaths_plot_std <- age_prop_deaths_plotdat %>%
  dplyr::filter(care_home_deaths=="yes") %>%
  ggplot() + theme_bw() +
  geom_line(aes(x = age_mid, y = cum_prop_deaths_std, group=study_id,color=study_id), size = 0.3) +
  scale_color_manual(values = col_vec, name = "study_id") +
  xlab("Age (yrs)") + ylab("Cumulative proportion of deaths, age-standardised") +
  xyaxis_plot_theme
if(write2file) ggsave(filename = "results/descriptive_figures/age_prop_deaths_plot_std.tiff", plot = age_prop_deaths_plot_std, width = 7, height = 5)




#......................
# rgn adj seroprevalence
#......................
rgn_seroplot <- rgnplotdat %>%
  dplyr::select(c("study_id", "region", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_point(aes(x = region, y = seroprevadj, color = study_id), size = 2.5) +
  scale_color_manual("Study ID", values = discrete_colors) +
  facet_wrap(.~study_id, scales = "free_x") +
  xlab("Region") + ylab("Adj. Seroprevalence (%)") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/rgn_adj_seroplot.jpg",
                           plot = rgn_seroplot)
#......................
# crude raw IFR
#......................
rgn_IFR_plot <- rgnplotdat %>%
  dplyr::filter(seromidpt == obsday) %>%
  dplyr::select(c("study_id", "region", "cumdeaths", "popn", "seroprev")) %>%
  dplyr::mutate(infxns = popn * seroprev,
                crudeIFR =  cumdeaths/(infxns+cumdeaths),
                crudeIFR = ifelse(crudeIFR > 1, 1, crudeIFR),
                seroprev = seroprev * 100 ) %>%
  dplyr::filter(infxns > 0) %>%
  ggplot() +
  geom_point(aes(x = region, y = crudeIFR, color = seroprev), size = 2.5) +
  facet_wrap(.~study_id, scales = "free_x") +
  scale_color_gradientn("Raw Seroprevalence (%)",
                        colors = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))) +
  xlab("Region") + ylab("Crude Infection Fatality Rate") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/rgn_IFR_raw_plot.jpg",
                           plot = rgn_IFR_plot)

#......................
# crude adj IFR
#......................
rgn_IFR_plot <- rgnplotdat %>%
  dplyr::filter(seromidpt == obsday) %>%
  dplyr::select(c("study_id", "region", "cumdeaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(infxns = popn * seroprevadj,
                crudeIFR =  cumdeaths/(infxns+cumdeaths),
                crudeIFR = ifelse(crudeIFR > 1, 1, crudeIFR),
                seroprevadj = seroprevadj * 100 ) %>%
  dplyr::filter(infxns > 0) %>%
  ggplot() +
  geom_point(aes(x = region, y = crudeIFR, color = seroprevadj), size = 2.5) +
  facet_wrap(.~study_id, scales = "free_x") +
  scale_color_gradientn("Adj. Seroprevalence (%)",
                        colors = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))) +
  xlab("Region") + ylab("Crude Infection Fatality Rate") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/rgn_IFR_adj_plot.jpg",
                           plot = rgn_IFR_plot)

#......................
# standardized deaths by seroprev
#......................
std_deaths_seroplotdat <- rgnplotdat %>%
  dplyr::filter(seromidpt == obsday)
# write out for later use
write.csv(std_deaths_seroplotdat, file = "data/derived/region_summ_IFR.csv")

# standardized deaths
std_deaths_seroplot <- std_deaths_seroplotdat %>%
  dplyr::select(c("study_id", "region", "std_cum_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = seroprevadj, y = std_cum_deaths, fill = study_id), shape = 21, size = 2.5, stroke = 0.2) +
  scale_fill_manual(values = col_vec, name = "study_id") +
  xlab("Adjusted Seroprevalence (%).") + ylab("Cumulative Deaths per Million") +
  #  labs(caption = "Cumulative deaths per million at midpoint of seroprevalence study") +
  xyaxis_plot_theme
if(write2file) ggsave(filename = "results/descriptive_figures/std_deaths_rgn_seroplot.tiff", plot = std_deaths_seroplot, width = 7, height = 5)

std_rgn_ifr_seroplot <- std_deaths_seroplotdat %>%
  dplyr::select(c("study_id", "region", "std_cum_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() + theme_bw() +
  geom_point(aes(x = seroprevadj, y = std_cum_deaths/(10000*seroprevadj), fill = study_id), shape = 21, size = 2.5, stroke = 0.2) +
  scale_fill_manual(values = col_vec, name = "study_id") +
  xlab("Adjusted Seroprevalence (%)") + ylab("IFR (%)") +
  xyaxis_plot_theme
if(write2file) ggsave(filename = "results/descriptive_figures/std_rgn_ifr_seroplot.tiff", plot = std_rgn_ifr_seroplot, width = 7, height = 5)

# standardized deaths with names of regions (busy plot for internal)
std_deaths_seroplot_busy <- std_deaths_seroplotdat %>%
  dplyr::select(c("study_id", "region", "std_cum_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_point(aes(x = seroprevadj, y = std_cum_deaths, color = study_id), size = 2) +
  ggrepel::geom_text_repel(aes(x = seroprevadj, y = std_cum_deaths, label = region)) +
  scale_color_manual("Study ID", values = discrete_colors) +
  xlab("Adjusted Seroprevalence (%).") + ylab("Cumulative Deaths per Million") +
  xyaxis_plot_theme
if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/std_deaths_rgn_seroplot_busy.jpg",
                           plot = std_deaths_seroplot_busy, width_wide = 8, height_wide = 5.5)



std_deaths_seroplot <- rgnplotdat %>%
  dplyr::filter(seromidpt == obsday) %>%
  dplyr::select(c("study_id", "region", "std_cum_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_point(aes(x = seroprevadj, y = std_cum_deaths, color = study_id), size = 1.2) +
  ggrepel::geom_text_repel(aes(x = seroprevadj, y = std_cum_deaths, label = region), size = 2.5) +
  facet_wrap(.~study_id) +
  scale_color_manual("Study ID", values = discrete_colors) +
  xlab("Adj. Seroprevalence (%).") + ylab("Cum. Deaths per Million") +
  labs(caption = "Cumulative Deaths per Million at midpoint of Seroprevalence Study") +
  xyaxis_plot_theme
if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/std_deaths_seroplot_labeled.jpg",
                           plot = std_deaths_seroplot)


#......................
# standardized deaths by rgn
#......................
rgn_std_cum_deaths_plot <- rgnplotdat %>%
  dplyr::filter(seromidpt == obsday) %>%
  dplyr::select(c("study_id", "region", "std_cum_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_point(aes(x = region, y = std_cum_deaths, fill = seroprevadj), color = "#000000", size = 2.5, shape = 21, alpha = 0.8) +
  facet_wrap(.~study_id, scales = "free_x") +
  scale_color_manual("Study ID", values = discrete_colors) +
  scale_fill_gradientn("Adj. Seroprevalence (%)",
                       colors = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))) +
  xlab("Region") + ylab("Cum. Deaths per Million") +
  labs(caption = "Cumulative Deaths per Million at midpoint of Seroprevalence Study") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/rgn_std_cum_deaths_plot.jpg",
                           plot = rgn_std_cum_deaths_plot)
#......................
# daily standardized deaths by rgn
#......................
rgn_std_daily_deaths_plot <- rgnplotdat %>%
  dplyr::select(c("study_id", "obsday", "region", "region", "std_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_line(aes(x = obsday, y = std_deaths, color = region), alpha = 0.8, size = 1.2) +
  facet_wrap(.~study_id) +
  xlab("Obs. Day") + ylab("Daily Deaths per Million") +
  xyaxis_plot_theme +
  theme(legend.position = "none")
if(write2file) jpgsnapshot(outpath = "figures/descriptive_figures/rgn_std_daily_deaths_plot.jpg",
                           plot = rgn_std_daily_deaths_plot)

#......................
# population structure
#......................
populationdf <- readr::read_tsv("data/raw/population.tsv") %>%
  dplyr::select(-c("reference")) %>%
  dplyr::filter(age_breakdown==1 & !is.na(study_id) & study_id!="IRN1" & study_id!="KEN1") %>%
  dplyr::arrange(study_id,age_low,age_high) %>%
  dplyr::group_by(study_id,age_high) %>%
  dplyr::summarise(pop=sum(population)) %>%
  dplyr::ungroup()

pop_tot<-populationdf %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(tot_pop=sum(pop)) %>%
  dplyr::select(study_id,tot_pop) %>%
  dplyr::ungroup()

populationdf <- full_join(populationdf,pop_tot,by="study_id") %>%
  dplyr::mutate(prop_pop=pop/tot_pop,
                age_high=replace(age_high,age_high==999,100))

cumu <- populationdf %>%
  dplyr::group_by(study_id) %>%
  dplyr::arrange(study_id,age_high) %>%
  dplyr::summarise(cum_prop_pop=cumsum(prop_pop))
populationdf$cum_prop_pop<-cumu$cum_prop_pop

pop_age_plot <-ggplot(populationdf, aes(x = age_high, y = cum_prop_pop,group=study_id)) +
  geom_line(aes(color = study_id), alpha = 0.8, size = 1) +
  xlab("Age") + ylab("Cumulative proportion of population") +
  coord_cartesian(xlim = c(50,100), ylim=c(0.5,1))  +
  xyaxis_plot_theme #+
if(write2file) jpgsnapshot(outpath = "results/descriptive_figures/pop_cum_age_plot.jpg",
                           plot = pop_age_plot,width_wide = 8,height_wide = 5.5)

over80<-populationdf %>%
  dplyr::mutate(ageband=cut(age_high,breaks=c(-1,81,1000))) %>%
  dplyr::group_by(study_id,ageband) %>%
  dplyr::summarise(prop_pop=sum(prop_pop)) %>%
  dplyr::filter(ageband=="(81,1e+03]")

###### OVERALL IFRS
## TODO - CHECK CHE1-2 (SEEM LOW)
age_IFRraw_plot0<-read.csv("data/derived/age_summ_IFR.csv")
ifr0<-age_IFRraw_plot0 %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(n_deaths=sum(cumdeaths),
                   infxns=sum(infxns),
                   pop=sum(popn)) %>%
  dplyr::mutate(ifr=n_deaths/(infxns+n_deaths))
ifr0<-full_join(ifr0,over80,by="study_id")
ifr0$prop_pop[which(ifr0$study_id=="GBR3")]<-0.04454408
ifr0$prop_pop[which(ifr0$study_id=="BRA1")]<-0.021
ifr0<-full_join(ifr0,prop_deaths_70,by="study_id")


par(mfrow=c(1,1))
plot(ifr0$prop_pop,ifr0$ifr*100,ylab="IFR (%)",xlab="proportion of population over 80",pch=19)


#............................................................
#---- Figure of Seroprevalence and Seroreversion #----
#...........................................................
datmap <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS")

# SeroPrevalences by age portion
SeroPrevPlotDat <- datmap %>%
  dplyr::filter(breakdown == "ageband") %>%
  dplyr::filter(!grepl("_nch", study_id)) %>%
  dplyr::select(c("study_id", "seroprev_adjdat")) %>%
  dplyr::filter(! study_id %in% c(c("CHE2", "DNK1", "LUX1", "NLD1",
                                  "SWE1", "LA_CA1"))) %>% # excluding studies w/ constant assumption
  tidyr::unnest(cols = "seroprev_adjdat")

# filter to latest date if multiple serosurveys
SeroPrevPlotDat <- SeroPrevPlotDat %>%
  dplyr::group_by(study_id, ageband) %>%
  dplyr::filter(obsdaymax == max(obsdaymax))

# add uncertainty in raw seroprevalence based on binomial
SeroPrevPlotDat_sub <- SeroPrevPlotDat %>%
  dplyr::filter(!is.na(n_positive)) %>%
  dplyr::filter(!is.na(n_tested)) %>%
  dplyr::mutate(crude_seroprev_obj = purrr::map2(n_positive, n_tested, .f = function(x,n){ binom.test(x,n) }),
                crude_seroprev_CI = purrr::map(crude_seroprev_obj, "conf.int"),
                crude_seroprevLCI = purrr::map_dbl(crude_seroprev_CI, function(x){x[[1]]}),
                crude_seroprevUCI = purrr::map_dbl(crude_seroprev_CI, function(x){x[[2]]}),
                crude_seroprev = purrr::map_dbl(crude_seroprev_obj, "estimate"))
# add back in ITA, which only have 95% CIs
SeroPrevPlotDat <- SeroPrevPlotDat %>%
  dplyr::filter(study_id == "ITA1") %>%
  dplyr::mutate(crude_seroprev_obj = NA,
                crude_seroprev_CI = NA,
                crude_seroprevLCI = serolci,
                crude_seroprevUCI = serouci,
                crude_seroprev = seroprev) %>%
  dplyr::bind_rows(., SeroPrevPlotDat_sub)

# plot out part A
FigA <- SeroPrevPlotDat %>%
  dplyr::mutate(age_low = as.numeric(stringr::str_split_fixed(ageband, "-", n = 2)[,1]),
                age_high = as.numeric(stringr::str_split_fixed(ageband, "-", n = 2)[,2]),
                age_high = ifelse(age_high == 999, 100, age_high),
                age_mid = (age_high + age_low)/2) %>%
  ggplot() +
  geom_pointrange(aes(x = age_mid, y = crude_seroprev, ymin = crude_seroprevLCI, ymax = crude_seroprevUCI,
                      color = study_id), alpha = 0.8) +
  geom_line(aes(x = age_mid, y = seroprev, color = study_id),
            alpha = 0.8, size = 1.2, show.legend = F) +
  scale_color_manual("Study ID", values = mycolors) +
  xlab("Age (yrs).") + ylab("Seroprevalence (%)") +
  xyaxis_plot_theme +
  theme(legend.position = "bottom") +
  theme(plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))

#............................................................
# SeroReversion Portion
#...........................................................
sero_rev_comb <- readRDS("results/sero_reversion/sero_rev_dat.RDS")
weibull_params <- readRDS("results/sero_reversion/weibull_params.RDS")
KM1_mod <- readRDS("results/sero_reversion/KaplanMeierFit.RDS")
survobj_km <- readRDS("results/sero_reversion/survobj_km.RDS")
WBmod <- readRDS("results/sero_reversion/WeibullFit.RDS")
serotime <- readRDS("results/sero_reversion/sero_reversion_incld_data.RDS")

#......................
# portion B Kaplan Meier
#......................
KMplot <- survminer::ggsurvplot(
  fit = KM1_mod,
  data = sero_rev_comb,
  ylab = "Prob. of Seropositivity Persistence",
  xlab = "Time (Days)",
  font.x = c(size = 12, face = "bold"),
  font.y = c(size = 12, face = "bold"),
  size = 1.5,
  conf.int = TRUE,
  conf.int.fill = "#bdbdbd",
  conf.int.alpha = 0.4,
  palette = "#3182bd", # color of line
  censor = TRUE, # add censoring points
  censor.shape = 124,
  censor.size = 5,
  ylim = c(0.5, 1))
FigB <- KMplot$plot +
  xyaxis_plot_theme +
  theme(axis.title = element_text(size = 12, face = "bold"),
        legend.position = "none",
        plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))



#......................
# portion C Weibull
#......................
# fitted 'survival'
t <- seq(from = 0, to = max(serotime$days_post_symptoms), by = 0.5)
cum_weib <- exp(-(t/weibull_params$wscale)^weibull_params$wshape)   # cumulative weibull
# mu is approx 129.62; sigma is 39.31, so 2*mean
d_weib <- dweibull(x = seq(0, 260, by = 0.5),
                   shape = weibull_params$wshape, scale = weibull_params$wscale)
d_weib_df <- tibble::tibble(time = 1:length(d_weib)/2,
                            d_weib = d_weib)

FigC <- ggplot() +
  geom_histogram(data = sero_rev_comb,
                 aes(x = time_to_event, y = ..density.., fill = factor(status)),
                 position = "identity",
                 alpha = 0.5) +
  geom_line(data = d_weib_df,
            aes(x = time, y = d_weib),
            size = 1.2, alpha = 0.9, color = "#bdbdbd") +
  scale_fill_manual(values = c("#3BABFD", "#EA4335")) +
  ylab("Density") +
  xlab("Seroreversion Times (Days)") +
  xyaxis_plot_theme +
  theme(legend.position = "none",
        plot.margin = unit(c(0.05, 0.05, 0.05, 1),"cm"))


# bring together
rght <- cowplot::plot_grid(FigB, FigC,
                           ncol = 1, nrow = 2, align = "v",
                           labels = c("(B)", "(C)"))

mainFig <- cowplot::plot_grid(FigA, rght,
                              ncol = 2, nrow = 1,
                              labels = c("(A)", ""))

jpeg("figures/final_figures/Figure_dscdat.jpg",
     width = 11, height = 8, units = "in", res = 500)
plot(mainFig)
graphics.off()





#............................................................
#---- Figure of Seroprevalence vs. Not-Modelled Adj. IFR #----
#...........................................................
source("R/delta_method.R")
datmap <- readRDS("results/descriptive_results/descriptive_results_datamap.RDS")

# get regions
rgns <- datmap %>%
  dplyr::filter(breakdown == "region") %>%
  dplyr::select(-c("care_home_deaths", "data", "seroprev_adjdat", "std_deaths")) %>%
  tidyr::unnest(cols = "plotdat")

#......................
# get standard errors
#......................
# Delta method needs standard error of seroprev SE(p)
# where SE(p) is the standard error of the binomial proportion for all studies except ITA
# for ITA, we use SE(p) as the logit transformed SE gleaned from the provided CIs
#
rgnsSE_sub <- rgns %>%
  dplyr::filter(study_id != "ITA1") %>%
  dplyr::group_by(study_id) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::ungroup(.) %>%
  dplyr::select(c("study_id", "region", "n_positive", "n_tested")) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::mutate(seroprev = n_positive/n_tested,
                binom_se = sqrt(seroprev * (1-seroprev))) %>%
  dplyr::select(c("study_id", "region", "binom_se"))

rgnsSE_ita <- rgns %>%
  dplyr::filter(study_id == "ITA1") %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::select(c("study_id", "region", "serolci", "serouci")) %>%
  dplyr::filter(!duplicated(.)) %>%
  dplyr::mutate(binom_se = (COVIDCurve:::logit(serouci) - COVIDCurve:::logit(serolci))/(1.96 * 2))  %>%
  dplyr::select(c("study_id", "region", "binom_se"))

rgnsSE <- dplyr::bind_rows(rgnsSE_sub, rgnsSE_ita)


#......................
# subset regions to parts we need
# and perform calculation
#......................
delta_IFR <- rgns %>%
  dplyr::group_by(study_id, region) %>%
  dplyr::filter(seromidpt == max(seromidpt)) %>% # latest serostudy
  dplyr::filter(obsday == seromidpt) %>%  # sero obs day
  dplyr::select(c("study_id", "region", "cumdeaths", "popn", "seroprev")) %>%
  dplyr::left_join(., rgnsSE, by = c("study_id", "region")) %>%
  dplyr::mutate(IFRcalc = cumdeaths  / (seroprev * popn + cumdeaths),
                crit = purrr::map_dbl(seroprev, get_delta_CI_val, deaths = cumdeaths, popN = popn, SE = binom_se),
                IFRbound = purrr::map(crit, getCI_from_logit_transfrom, pt = IFRcalc, alpha = 1.96, tol = 1e-3),
                lower_ci = purrr::map_dbl(IFRbound, "lower.ci"),
                upper_ci = purrr::map_dbl(IFRbound, "upper.ci"))

#......................
# make plots
#......................
