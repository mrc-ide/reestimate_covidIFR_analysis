#................................................................................................
## Purpose: Plot descriptive statistics
##
## Notes:
#................................................................................................
#......................
# discrete colors inspired by https://www.schemecolor.com/
#......................
discrete_colors <- c("#CB00FF", # vivid orchid
                     "#E59500", # harvest gold
                     "#002642", # oxford blue
                     "#761EF8", # electric indigo
                     "#3EF944", # neon green
                     "#F61843", # crayola red
                     "#271F8E", # Cosmic Cobalt
                     "#5945AA", # plump purple
                     "#82AA96", # morning blue-green
                     "#66678E", # Dark Blue-Gray
                     "#840032", # burgundy
                     "#895B51" # spicy mix brown
                     )


#......................
# setup
#......................
library(tidyverse)
source("R/crude_plot_summ.R")
source("R/my_themes.R")
source("R/extra_plotting_functions.R")
dir.create("figures/descriptive_figures/", recursive = TRUE)

#............................................................
# read in data map
#...........................................................
datmap <- readxl::read_excel("data/derived/derived_data_map.xlsx")
datmap <- datmap %>%
  dplyr::mutate(data = purrr::map(relpath, readRDS))

## assign constant colours to each study
study_cols<-data.frame(study_id=unique(datmap$study_id),cols=discrete_colors[1:length(unique(datmap$study_id))])

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
    deathdat_long = purrr::map(data, "deathsMCMC"),
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
saveRDS(datmap, file = "data/derived/descriptive_results_datamap.RDS")


#............................................................
# Age Bands Plots/Descriptions
#...........................................................
ageplotdat <- datmap %>%
  dplyr::filter(breakdown == "ageband") %>%
  dplyr::select(c("study_id", "plotdat")) %>%
  tidyr::unnest(cols = "plotdat")
ageplotdat<-full_join(ageplotdat,study_cols,by="study_id")

###### filter to only plot the latest serology when there are multiple rounds - ok?
maxDays <- ageplotdat %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(max_day=max(obsdaymax))
ageplotdat<-full_join(ageplotdat,maxDays,by="study_id")
ageplotdat<-filter(ageplotdat,obsdaymax==max_day)


#......................
# age adj seroprevalence
#......................
age_seroplot <- ageplotdat %>%
  dplyr::select(c("study_id", "age_mid", "seroprev")) %>%
  dplyr::mutate(seroprev = seroprev * 100) %>%
  ggplot() +
  geom_line(aes(x = age_mid, y = seroprev, color = study_id), alpha = 0.8, size = 1.2) +
  geom_point(aes(x = age_mid, y = seroprev, color = study_id)) +
  scale_color_manual("Study ID", values = discrete_colors) +
  xlab("Age (yrs).") + ylab("Raw Seroprevalence (%)") +
  xyaxis_plot_theme
jpgsnapshot(outpath = "figures/descriptive_figures/age_raw_seroplot.jpg",
            plot = age_seroplot)

#......................
# age adj seroprevalence
#......................
age_seroplot <- ageplotdat %>%
  dplyr::select(c("study_id", "age_mid", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_line(aes(x = age_mid, y = seroprevadj, color = study_id), alpha = 0.8, size = 1.2) +
  geom_point(aes(x = age_mid, y = seroprevadj, color = study_id)) +
  scale_color_manual("Study ID", values = discrete_colors) +
  xlab("Age (yrs).") + ylab("Adj. Seroprevalence (%)") +
  xyaxis_plot_theme
jpgsnapshot(outpath = "figures/descriptive_figures/age_adj_seroplot.jpg",
            plot = age_seroplot)

#......................
# crude IFR
#......................
# raw serology
age_IFRraw_plot0 <- ageplotdat %>%
  dplyr::filter(seromidpt == obsday) %>%
  dplyr::select(c("study_id","n_positive","n_tested","ageband","age_low","age_high", "age_mid", "cumdeaths", "popn", "seroprev", "seroprevadj")) %>%
  dplyr::mutate(infxns = popn * seroprev,
                crudeIFR =  cumdeaths/(infxns+cumdeaths),
                crudeIFR = ifelse(crudeIFR > 1, 1, crudeIFR),
                seroprev = seroprev * 100 )

# write this out for later use
readr::write_csv(age_IFRraw_plot0, path = "data/derived/age_summ_IFR.csv")

age_IFRraw_plot <- age_IFRraw_plot0 %>%
  ggplot() +
  geom_line(aes(x = age_mid, y = crudeIFR, color = study_id), alpha = 0.8, size = 1.2) +
  geom_point(aes(x = age_mid, y = crudeIFR, fill = seroprev), color = "#000000", size = 2.5, shape = 21, alpha = 0.8) +
  scale_color_manual("Study ID", values = discrete_colors) +
  scale_fill_gradientn("Raw Seroprev.",
                       colors = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))) +
  xlab("Age (yrs).") + ylab("Crude Infection Fatality Rate") +
  xyaxis_plot_theme
jpgsnapshot(outpath = "figures/descriptive_figures/age_IFRraw_plot.jpg",
            plot = age_IFRraw_plot)

### probably a neater way to do this. TODO
age_IFRraw_plot2 <- age_IFRraw_plot0 %>%
  ggplot() +
  geom_point(aes(x = age_mid, y = crudeIFR, color = study_id))+ #, color = "#000000", size = 2.5, shape = 21, alpha = 0.8) +
  geom_line(aes(x = age_mid, y = crudeIFR, color = study_id), alpha = 0.8, size = 1.2) +
  scale_color_manual("Study ID", values = discrete_colors) +
  xlab("Age (years)") + ylab("Crude infection fatality rate") +
  xyaxis_plot_theme
jpgsnapshot(outpath = "figures/descriptive_figures/age_IFRraw_plot2.jpg",
            plot = age_IFRraw_plot2, width_wide = 8,height_wide = 5.5)

age_IFRraw_plot_log <- age_IFRraw_plot0 %>%
  filter(crudeIFR>0) %>%
  ggplot() +
  geom_point(aes(x = age_mid, y = crudeIFR, color = study_id))+ #, color = "#000000", size = 2.5, shape = 21, alpha = 0.8) +
  geom_line(aes(x = age_mid, y = crudeIFR, color = study_id), alpha = 0.8, size = 1.2) +
  scale_color_manual("Study ID", values=discrete_colors) +
  xlab("Age (years)") + ylab("Crude infection fatality rate") +
  xyaxis_plot_theme +
  scale_y_log10() +
#  coord_cartesian(ylim=c(0.00000000001,1))
jpgsnapshot(outpath = "results/descriptive_figures/age_IFRraw_plot_log.jpg",
            plot = age_IFRraw_plot_log,width_wide = 8,height_wide = 5.5)

# seroadj
age_IFRadj_plot <- ageplotdat %>%
  dplyr::filter(seromidpt == obsday) %>%
  dplyr::select(c("study_id", "age_mid", "cumdeaths", "popn", "seroprev", "seroprevadj")) %>%
  dplyr::mutate(infxns = popn * seroprevadj,
                crudeIFR =  cumdeaths/(cumdeaths+infxns),
                crudeIFR = ifelse(crudeIFR > 1, 1, crudeIFR),
                seroprevadj = seroprevadj * 100 )

%>%
  ggplot() +
  geom_line(aes(x = age_mid, y = crudeIFR, color = study_id), alpha = 0.8, size = 1.2) +
  geom_point(aes(x = age_mid, y = crudeIFR, fill = seroprevadj), color = "#000000", size = 2.5, shape = 21, alpha = 0.8) +
  scale_color_manual("Study ID", values = discrete_colors) +
  scale_fill_gradientn("Adj. Seroprev.",
                       colors = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))) +
  xlab("Age (yrs).") + ylab("Crude Infection Fatality Rate") +
  xyaxis_plot_theme
jpgsnapshot(outpath = "figures/descriptive_figures/age_IFRadj_plot.jpg",
            plot = age_IFRadj_plot)

#......................
# standardized deaths by seroprev
#......................
std_deaths_seroplot <- ageplotdat %>%
  dplyr::filter(seromidpt == obsday) %>%
  dplyr::select(c("study_id", "age_mid", "std_cum_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_point(aes(x = seroprevadj, y = std_cum_deaths, color = study_id, size = age_mid), size = 1.2) +
  scale_color_manual("Study ID", values = discrete_colors) +
  scale_size("Mid. Age", range = c(0.1, 3)) +
  xlab("Adj. Seroprevalence (%).") + ylab("Cum. Deaths per Million") +
  labs(caption = "Cumulative Deaths per Million at midpoint of Seroprevalence Study") +
  xyaxis_plot_theme
jpgsnapshot(outpath = "figures/descriptive_figures/std_deaths_age_seroplot.jpg", # had same name as regional data plot previously
            plot = std_deaths_seroplot)

#......................
# standardized deaths by age
#......................
age_std_cum_deaths_plot <- ageplotdat %>%
  dplyr::filter(seromidpt == obsday) %>%
  dplyr::select(c("study_id", "age_mid", "std_cum_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_line(aes(x = age_mid, y = std_cum_deaths, color = study_id), alpha = 0.8, size = 1.2) +
  geom_point(aes(x = age_mid, y = std_cum_deaths, fill = seroprevadj), color = "#000000", size = 2.5, shape = 21, alpha = 0.8) +
  scale_color_manual("Study ID", values = discrete_colors) +
  scale_fill_gradientn("Adj. Seroprevalence (%)",
                       colors = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous"))) +
  xlab("Age (yrs).") + ylab("Cum. Deaths per Million") +
  labs(caption = "Cumulative Deaths per Million at midpoint of Seroprevalence Study") +
  xyaxis_plot_theme
jpgsnapshot(outpath = "figures/descriptive_figures/age_std_cum_deaths_plot.jpg",
            plot = age_std_cum_deaths_plot)


#......................
# cumulative proportion deaths by age
#......................
age_IFRraw_plot0 <- age_IFRraw_plot0 %>%
  dplyr::mutate(d_per_mill=cumdeaths/popn)

tot_deaths <- age_IFRraw_plot0 %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(tot_deaths=sum(cumdeaths),
                   tot_deaths_std=sum(d_per_mill)) %>%
  dplyr::select(study_id, tot_deaths,tot_deaths_std) %>%
  ungroup()

age_prop_deaths_plot<- full_join(age_IFRraw_plot0,tot_deaths,by="study_id") %>%
  dplyr::mutate(prop_deaths = cumdeaths/tot_deaths,
                prop_deaths_std = d_per_mill/tot_deaths_std) %>%
  dplyr::arrange(study_id,age_mid)

cumu_deaths<- age_prop_deaths_plot %>%
  dplyr::group_by(study_id) %>%
  dplyr::summarise(cum_prop_deaths=cumsum(prop_deaths),
                   cum_prop_deaths_std=cumsum(prop_deaths_std))

age_prop_deaths_plot$cum_prop_deaths<-cumu_deaths$cum_prop_deaths
age_prop_deaths_plot$cum_prop_deaths_std<-cumu_deaths$cum_prop_deaths_std


########## Raw deaths by age cumulative (just showing us the population structure more than anything?)
age_prop_deaths_plot<-ggplot(age_prop_deaths_plot, aes(x = age_mid, y = cum_prop_deaths, group=study_id)) +
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
jpgsnapshot(outpath = "results/descriptive_figures/age_prop_cum_deaths_plot.jpg",
            plot = age_prop_deaths_plot)

########## Deaths per capita by age cumulative
age_prop_deaths_plot_std<-ggplot(age_prop_deaths_plot, aes(x = age_mid, y = cum_prop_deaths_std, group=study_id)) +
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
  xlab("Age (yrs)") + ylab("Cumulative proportion of deaths, age-standardised") +
  xyaxis_plot_theme
jpgsnapshot(outpath = "results/descriptive_figures/age_prop_cum_deaths_std_plot.jpg",
            plot = age_prop_deaths_plot_std,width_wide = 8,height_wide = 5.5)


#......................
# daily standardized deaths by age
#......................
age_std_daily_deaths_plot <- ageplotdat %>%
  dplyr::select(c("study_id", "obsday", "ageband", "age_mid", "std_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(ageband = forcats::fct_reorder(ageband, age_mid),
                seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_line(aes(x = obsday, y = std_deaths, color = ageband), alpha = 0.8, size = 1.2) +
  facet_wrap(.~study_id, scales = "free_y") +
  xlab("Obs. Day") + ylab("Daily Deaths per Million") +
  xyaxis_plot_theme
jpgsnapshot(outpath = "figures/descriptive_figures/age_std_daily_deaths_plot.jpg",
            plot = age_std_daily_deaths_plot)


#............................................................
# Regional Plots/Descriptions
#............................................................
rgnplotdat <- datmap %>%
  dplyr::filter(breakdown == "region") %>%
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
rgn_seroplot <- rgnplotdat %>%
  dplyr::select(c("study_id", "region", "seroprev")) %>%
  dplyr::mutate(seroprev = seroprev * 100) %>%
  ggplot() +
  geom_point(aes(x = region, y = seroprev, color = study_id), size = 2.5) +
  scale_color_manual("Study ID", values = discrete_colors) +
  facet_wrap(.~study_id, scales = "free_x") +
  xlab("Region") + ylab("Raw Seroprevalence (%)") +
  xyaxis_plot_theme +
  theme(axis.text.x = element_text(family = "Helvetica", hjust = 1, size = 8, angle = 45))

jpgsnapshot(outpath = "figures/descriptive_figures/rgn_raw_seroplot.jpg",
            plot = rgn_seroplot)



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

jpgsnapshot(outpath = "figures/descriptive_figures/rgn_adj_seroplot.jpg",
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

jpgsnapshot(outpath = "figures/descriptive_figures/rgn_IFR_raw_plot.jpg",
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

jpgsnapshot(outpath = "figures/descriptive_figures/rgn_IFR_adj_plot.jpg",
            plot = rgn_IFR_plot)

#......................
# standardized deaths by seroprev
#......................
std_deaths_seroplot <- rgnplotdat %>%
  dplyr::filter(seromidpt == obsday)
# write out for later use
readr::write_csv(std_deaths_seroplot, path = "data/derived/region_summ_IFR.csv")

std_deaths_seroplot <- std_deaths_seroplot %>%
  dplyr::select(c("study_id", "region", "std_cum_deaths", "popn", "seroprevadj")) %>%
  dplyr::mutate(seroprevadj = seroprevadj * 100) %>%
  ggplot() +
  geom_point(aes(x = seroprevadj, y = std_cum_deaths, color = study_id), size = 2) +
  scale_color_manual("Study ID", values = discrete_colors) +
  xlab("Adjusted Seroprevalence (%).") + ylab("Cumulative Deaths per Million") +
#  labs(caption = "Cumulative deaths per million at midpoint of seroprevalence study") +
  xyaxis_plot_theme
jpgsnapshot(outpath = "figures/descriptive_figures/std_deaths_rgn_seroplot.jpg",
            plot = std_deaths_seroplot,width_wide = 8,height_wide = 5.5)



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
jpgsnapshot(outpath = "figures/descriptive_figures/std_deaths_seroplot_labeled.jpg",
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

jpgsnapshot(outpath = "figures/descriptive_figures/rgn_std_cum_deaths_plot.jpg",
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
jpgsnapshot(outpath = "figures/descriptive_figures/rgn_std_daily_deaths_plot.jpg",
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
jpgsnapshot(outpath = "results/descriptive_figures/pop_cum_age_plot.jpg",
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
ifr0$prop_pop.y[which(ifr0$study_id=="GBR3")]<-0.04454408
ifr0$prop_pop.y[which(ifr0$study_id=="BRA1")]<-0.021

par(mfrow=c(1,1))
plot(ifr0$prop_pop.y,ifr0$ifr*100,ylab="IFR (%)",xlab="proportion of population over 80",pch=19)
