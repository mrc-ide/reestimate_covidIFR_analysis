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
library(pals)

#............................................................
# Age Bands
#...........................................................
# Spain
esp_orig <- get_orig_seroprev(IFRmodinput = readRDS("data/derived/ESP/ESP_agebands.RDS"),
                              groupingvar = "ageband") %>%
  dplyr::mutate(country = "ESP")

# Netherlands
nld_orig <- get_orig_seroprev(IFRmodinput = readRDS("data/derived/NLD/NLD_agebands.RDS"),
                              groupingvar = "ageband") %>%
  dplyr::mutate(country = "NLD")
# Denmark
dnk_orig <- get_orig_seroprev(IFRmodinput = readRDS("data/derived/DNK/DNK_agebands.RDS"),
                              groupingvar = "ageband") %>%
  dplyr::mutate(country = "DNK")
# Switzerland
che_orig <- get_orig_seroprev(IFRmodinput = readRDS("data/derived/CHE/CHE_agebands.RDS"),
                              groupingvar = "ageband") %>%
  dplyr::mutate(country = "CHE")

# Sweden
swe_orig <- get_orig_seroprev(IFRmodinput = readRDS("data/derived/SWE/SWE_agebands.RDS"),
                              groupingvar = "ageband") %>%
  dplyr::mutate(country = "SWE")

#......................
# more in depth
#......................
# Spain
esp_res <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/ESP/ESP_agebands.RDS"),
                        groupingvar = "ageband") %>%
  dplyr::mutate(country = "ESP")

# Netherlands
nld_res <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/NLD/NLD_agebands.RDS"),
                        groupingvar = "ageband") %>%
  dplyr::mutate(country = "NLD")
# Denmark
dnk_res <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/DNK/DNK_agebands.RDS"),
                        groupingvar = "ageband") %>%
  dplyr::mutate(country = "DNK")
# Switzerland
che_res <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/CHE/CHE_agebands.RDS"),
                        groupingvar = "ageband") %>%
  dplyr::mutate(country = "CHE")

# Sweden
swe_res <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/SWE/SWE_agebands.RDS"),
                        groupingvar = "ageband") %>%
  dplyr::mutate(country = "SWE")

#............................................................
# Regional
#...........................................................
# Spain
esp_resr <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/ESP/ESP_regions.RDS"),
                         groupingvar = "region") %>%
  dplyr::mutate(country = "ESP")

# Netherlands
nld_resr <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/NLD/NLD_regions.RDS"),
                         groupingvar = "region") %>%
  dplyr::mutate(country = "NLD")

# Denmark
dnk_resr <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/DNK/DNK_regions.RDS"),
                         groupingvar = "region") %>%
  dplyr::mutate(country = "DNK")

# Sweden
swe_resr <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/SWE/SWE_regions.RDS"),
                         groupingvar = "region") %>%
  dplyr::mutate(country = "SWE")

# USA - Los Angeles, California
la_ca_resr <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/USA/LA_CA_regions.RDS"),
                           groupingvar = "region") %>%
  dplyr::mutate(country = "LA_CA")

# USA - Santa Clara, California
sc_ca_resr <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/USA/SC_CA_regions.RDS"),
                           groupingvar = "region") %>%
  dplyr::mutate(country = "SC_CA")

# USA - CH_MA
ch_ma_resr <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/USA/CH_MA_regions.RDS"),
                           groupingvar = "region") %>%
  dplyr::mutate(country = "CH_MA")

# USA - MD_FL
md_fl_resr <- get_crude_summarydf(IFRmodinput = readRDS("data/derived/USA/MD_FL_regions.RDS"),
                           groupingvar = "region") %>%
  dplyr::mutate(country = "MD_FL")

#............................................................
# plots
#...........................................................
#......................
# seroprevalence by age
#......................
agedat <- rbind.data.frame(esp_orig, nld_orig, dnk_orig, che_orig, swe_orig)

agedat %>%
  ggplot() +
  geom_line(aes(x = age_mid, y = seroprev, color = country), alpha = 0.8) +
  geom_point(aes(x = age_mid, y = seroprev, color = country)) +
  scale_color_manual(values = wesanderson::wes_palette("Darjeeling1")) +
  xlab("Age (yrs.)") + ylab("Seroprevalnce") +
  theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
        axis.text = element_text(family = "Helvetica", hjust = 0.5, size = 11),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(family = "Helvetica", hjust = 0.5, vjust = 0.5, size = 10),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "#000000", size = 1))





#......................
# DEATHS BY AGE
#......................
deaths<-read.csv("data/raw/deaths.csv")
deaths<- deaths %>%
  dplyr::filter(age_breakdown==1) %>%
  dplyr::mutate(age_high=recode(age_high,`999`=99L),
                age_mid=0.5*(age_low+age_high)) %>%
  dplyr::group_by(study_id,age_mid) %>%
  dplyr::summarise(n_deaths=sum(n_deaths)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(study_id) %>% ### add tot deaths per study.
  dplyr::mutate(deaths_denom=sum(n_deaths),
                deaths.prop.age=n_deaths/deaths_denom,
                cols=as.factor(study_id))
levels(deaths$cols)<-cols25(25)


par(mar=c(5,4,4,5))
plot(deaths$age_mid,deaths$deaths.prop.age,col=deaths$cols,xlab="age",ylab="proportion of deaths")
legend(102,0.7,unique(deaths$study_id),xpd=T,pch=1,col=unique(deaths$cols))

### Deaths per population by age
col_vec <- RColorBrewer::brewer.pal(7, "Set1")
names(col_vec) <- c("Spain", "Sweden", "Switzerland", "Denmark","Netherlands","United Kingdom","United States - LA")
par(mar=c(5,4,4,2))
scale<-1
plot(esp_res$age_mid,scale*esp_res$prop_deaths_per_pop,xlab="age group (years)",ylab="deaths per capita",xlim=c(0,100),
     pch=21, col.main="black", bg=col_vec[1],ylim=c(0,1))
points(che_res$age_mid,scale*che_res$prop_deaths_per_pop,pch=21, col.main="black", bg=col_vec[3])
points(dnk_res$age_mid,scale*dnk_res$prop_deaths_per_pop,pch=21, col.main="black", bg=col_vec[4])
points(nld_res$age_mid,scale*nld_res$prop_deaths_per_pop,pch=21, col.main="black", bg=col_vec[5])
legend(0,1.4,names(col_vec)[c(1,3,4,5,7)],pch=rep(21,4),col=rep("black",4), bty='n',
       pt.bg=col_vec[c(1,3,4,5,7)],xpd=T,ncol=2)



######## IFR BY AGE
col_vec <- RColorBrewer::brewer.pal(7, "Set1")
names(col_vec) <- c("Spain", "Sweden", "Switzerland", "Denmark","Netherlands","United Kingdom","United States - LA")
par(mar=c(5,4,4,2))
plot(esp_res$age_mid,100*esp_res$ifr_age_crude,xlab="age group (years)",ylab="Crude IFR",xlim=c(0,100),
     ylim=c(0,17),pch=21, col.main="black", bg=col_vec[1])
points(che_res$age_mid,100*che_res$ifr_age_crude,pch=21, col.main="black", bg=col_vec[3])
points(dnk_res$age_mid,100*dnk_res$ifr_age_crude,pch=21, col.main="black", bg=col_vec[4])
points(nld_res$age_mid,100*nld_res$ifr_age_crude,pch=21, col.main="black", bg=col_vec[5])
points(la_ca_res$age_mid,100*la_ca_res$ifr_age_crude,pch=21, col.main="black", bg=col_vec[7])
legend(0,23,names(col_vec)[c(1,3,4,5,7)],pch=rep(21,4),col=rep("black",4), bty='n',
       pt.bg=col_vec[c(1,3,4,5,7)],xpd=T,ncol=2)


######## REGION VERSUS MORTALITY
tiff(file="figures/sero_vs_deaths_region.tiff", width=2000,height=2200,res=300,compression="lzw")
col_vec <- c(RColorBrewer::brewer.pal(7, "Set1"),cols25(25))
legend_text <- c("Spain", "Sweden", "Switzerland", "Denmark","Netherlands","UK","US - LA",
                 "US - Santa Clara", "US - Massachusetts", "US - Miami-Dade", "US - NYC")
par(mar=c(5,4,9,2))
plot(esp_resr$seroprevalence,esp_resr$deaths_per_million,xlab="Seroprevalence (%)",ylab="deaths per million",
     ylim=c(0,1800),xlim=c(0,0.32),pch=21, col.main="black", bg=col_vec[1])
studies<-c("swe","che","dnk","nld","gbr","la_ca","sc_ca","ch_ma","md_fl","nyc_resr")
points(che_resr$seroprev,che_resr$deaths_per_million,pch=21, col.main="black", bg=col_vec[3])
points(dnk_resr$seroprevalence,dnk_resr$deaths_per_million,pch=21, col.main="black", bg=col_vec[4])
points(nld_resr$seroprevalence,nld_resr$deaths_per_million,pch=21, col.main="black", bg=col_vec[5])
points(la_ca_resr$seroprevalence,la_ca_resr$deaths_per_million,pch=21, col.main="black", bg=col_vec[11])
points(sc_ca_resr$seroprevalence,sc_ca_resr$deaths_per_million,pch=21, col.main="black", bg=col_vec[12])
points(ch_ma_resr$seroprevalence,ch_ma_resr$deaths_per_million,pch=21, col.main="black", bg=col_vec[13])
points(md_fl_resr$seroprevalence,md_fl_resr$deaths_per_million,pch=21, col.main="black", bg=col_vec[14])
points(nyc_resr$seroprevalence,nyc_resr$deaths_per_million,pch=21, col.main="black", bg=col_vec[15])
incl<-c(1,3,4,5,7:11)
legend(0,2400,legend_text[incl],pch=rep(21,4),col=rep("black",4), bty='n',
       pt.bg=col_vec[c(1,3,4,5,11:15)],xpd=T,ncol=2)
dev.off()
