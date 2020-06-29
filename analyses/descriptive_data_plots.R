########################
# Plot descriptive statistics
########################

library(dplyr)
library(pals)
#source("analyses/run_process_country_data.R")
rogan_gladen<-function(obs_prev,sens,spec) (obs_prev + spec -1)/(spec+sens-1)
sero_sheet<-read.csv("data/raw/seroprevalence.csv")

write2file<-F

#####################################################
## Align serology and deaths by age where possible. Document assumptions where no perfect alignment
#####################################################
### Spain. Perfect alignment of deaths and seroprevalence age groups.
i<-"ESP"
x<-readRDS(paste0("data/derived/",i,"/",i,"_agebands.RDS"))
curr_sero<-x$seroprev_group
curr_sero<-curr_sero %>%
  dplyr::mutate(ageband_deaths=cut(curr_sero$age_high, breaks=c(0,x$deaths_group$age_high))) %>%
  dplyr::group_by(ageband_deaths) %>%
  dplyr::summarise(ObsDaymin=mean(ObsDaymin),
                   ObsDaymax=mean(ObsDaymax),
                   n_tested=sum(n_tested),
                   n_positive=sum(n_positive)) %>%
  dplyr::mutate(seroprevalence=n_positive/n_tested) %>%
  dplyr::ungroup()

curr_sero$age_mid<-0.5*(x$deaths_group$age_low + x$deaths_group$age_high)
curr_sero$age_mid[nrow(curr_sero)]<-95
curr_sero$pop<-x$popN*x$prop_pop$pop_prop
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_crude<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$seroprev_adj_ss<-rogan_gladen(curr_sero$seroprevalence, x$sero_sens,x$sero_spec)
curr_sero$inf_pop_adj_ss<-curr_sero$seroprev_adj_ss * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_adj_ss<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_adj_ss
curr_sero$deaths_per_pop<-x$deaths_group$deaths_at_sero / curr_sero$pop
curr_sero$prop_deaths_per_pop<-curr_sero$deaths_per_pop/sum(curr_sero$deaths_per_pop)
esp_res<-curr_sero

### Netherlands. Not perfectly aligned.
# Assumptions.
# 1) 31-40 seroprevalence will be equivalent to the 30-39 age group, same for other 10 year age bands 30-70.
# 2) 18-30 seroprevalence = 20-29 seroprevalence
# 3) 60-72 seroprevalence = 60-69 seroprevalence
# 4) all other seroprevalence = national average.
i<-"NLD"
x<-readRDS(paste0("data/derived/",i,"/",i,"_agebands.RDS"))
curr_sero<-x$deaths_group
curr_sero$age_mid<-0.5*(curr_sero$age_low + curr_sero$age_high)
curr_sero$age_mid[which(curr_sero$ageband=='94-999')]<-98  ## TODO check what average is in this age group.
## enter national average by default.
curr_sero$seroprevalence<-x$seroprev$seroprev
#######
# enter for specific age groups
i<-which(x$seroprev_group$ageband=='18-30')
curr_sero$seroprevalence[which(curr_sero$ageband=='19-24' | curr_sero$ageband=='24-29')]<-
      x$seroprev_group$seroprevalence[i]
i<-which(x$seroprev_group$ageband=='30-40')
curr_sero$seroprevalence[which(curr_sero$ageband=='29-34' | curr_sero$ageband=='34-39')]<-
  x$seroprev_group$seroprevalence[i]
i<-which(x$seroprev_group$ageband=='40-50')
curr_sero$seroprevalence[which(curr_sero$ageband=='39-44' | curr_sero$ageband=='44-49')]<-
  x$seroprev_group$seroprevalence[i]
i<-which(x$seroprev_group$ageband=='50-60')
curr_sero$seroprevalence[which(curr_sero$ageband=='49-54' | curr_sero$ageband=='54-59')]<-
  x$seroprev_group$seroprevalence[i]
i<-which(x$seroprev_group$ageband=='60-72')
curr_sero$seroprevalence[which(curr_sero$ageband=='59-64' | curr_sero$ageband=='64-69')]<-
  x$seroprev_group$seroprevalence[i]

x$prop_pop<-x$prop_pop[which(!is.na(x$prop_pop$ageband)),]

curr_sero$pop<-x$popN*x$prop_pop$pop_prop
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_crude<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$seroprev_adj_ss<-rogan_gladen(curr_sero$seroprevalence, x$sero_sens,x$sero_spec)
curr_sero$inf_pop_adj_ss<-curr_sero$seroprev_adj_ss * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_adj_ss<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_adj_ss
curr_sero$deaths_per_pop<-x$deaths_group$deaths_at_sero / curr_sero$pop
curr_sero$prop_deaths_per_pop<-curr_sero$deaths_per_pop/sum(curr_sero$deaths_per_pop)
nld_res<-curr_sero


### Denmark
i<-"DNK"
x<-readRDS(paste0("data/derived/",i,"/",i,"_agebands.RDS"))
curr_sero<-x$deaths_group
curr_sero$age_mid<-0.5*(curr_sero$age_low + curr_sero$age_high)
curr_sero$age_mid[which(curr_sero$ageband=='89-999')]<-95
## enter national average by default.
curr_sero$seroprevalence<-x$seroprev$seroprev
curr_sero$seroprevalence[which(curr_sero$ageband=='59-69')]<-
  x$seroprev_group$seroprevalence[which(x$seroprev_group$ageband=='59-69')]

curr_sero$pop<-x$popN*x$prop_pop$pop_prop
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_crude<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$seroprev_adj_ss<-rogan_gladen(curr_sero$seroprevalence, x$sero_sens,x$sero_spec)
curr_sero$inf_pop_adj_ss<-curr_sero$seroprev_adj_ss * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_adj_ss<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_adj_ss
curr_sero$deaths_per_pop<-x$deaths_group$deaths_at_sero / curr_sero$pop
curr_sero$prop_deaths_per_pop<-curr_sero$deaths_per_pop/sum(curr_sero$deaths_per_pop)
dnk_res<-curr_sero

#################################
### Switzerland. Age bands not perfectly aligned
# Assume 5-19 seroprevalence = 0-20
i<-"CHE"
x<-readRDS(paste0("data/derived/",i,"/",i,"_agebands.RDS"))
curr_sero<-x$deaths_group
curr_sero$age_mid<-0.5*(curr_sero$age_low + curr_sero$age_high)
curr_sero$age_mid[which(curr_sero$ageband=='80-999')]<-90 ## TODO check exact.
## enter national average by default.
curr_sero$seroprevalence<-x$seroprev$seroprev
#### Enter for specific age groups
i<-which(x$seroprev_group$ageband=='5-19')
curr_sero$seroprevalence[which(curr_sero$ageband=='0-10' | curr_sero$ageband=='10-20')]<-
  x$seroprev_group$seroprevalence[i]
i<-which(x$seroprev_group$ageband=='19-49')
curr_sero$seroprevalence[which(curr_sero$ageband=='20-30' | curr_sero$ageband=='30-40'  |
                                 curr_sero$ageband=='40-50')]<-
    x$seroprev_group$seroprevalence[i]
i<-which(x$seroprev_group$ageband=='49-999')
curr_sero$seroprevalence[which(curr_sero$age_low>=50)]<-
  x$seroprev_group$seroprevalence[i]


curr_sero$pop<-x$popN*x$prop_pop$pop_prop
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_crude<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$seroprev_adj_ss<-rogan_gladen(curr_sero$seroprevalence, x$sero_sens,x$sero_spec)
curr_sero$inf_pop_adj_ss<-curr_sero$seroprev_adj_ss * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_adj_ss<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_adj_ss
curr_sero$deaths_per_pop<-x$deaths_group$deaths_at_sero / curr_sero$pop
curr_sero$prop_deaths_per_pop<-curr_sero$deaths_per_pop/sum(curr_sero$deaths_per_pop)
che_res<-curr_sero



#######################
#################################
### USA: Los Angeles county, California. Age bands not perfectly aligned
i<-"LA_CA"
x<-readRDS(paste0("data/derived/USA/",i,"_agebands.RDS"))
curr_sero<-x$deaths_group
curr_sero$age_mid<-0.5*(curr_sero$age_low + curr_sero$age_high)
curr_sero$age_mid[which(curr_sero$ageband=='65-999')]<-73 ## TODO check exact.
## enter average by default.
curr_sero$seroprevalence<-x$seroprev$seroprevalence

curr_sero$pop<-x$popN*x$prop_pop$pop_prop
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * x$popN*x$prop_pop$pop_prop
curr_sero$ifr_age_crude<-x$deaths_group$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$deaths_per_pop<-x$deaths_group$deaths_at_sero / curr_sero$pop
curr_sero$prop_deaths_per_pop<-curr_sero$deaths_per_pop/sum(curr_sero$deaths_per_pop)
la_ca_res<-curr_sero


#######################
# PROCESS REGIONAL INFO
### Spain.
i<-"ESP"
x<-readRDS(paste0("data/derived/",i,"/",i,"_regions.RDS"))
curr_sero<-x$deaths_group
curr_sero$seroprevalence<-x$seroprev_group$seroprevalence
curr_sero$pop<-x$popN*x$prop_pop$pop_prop
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * curr_sero$pop
curr_sero$ifr_crude<-curr_sero$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$deaths_per_million<-1000000*curr_sero$deaths_at_sero/curr_sero$pop
esp_resr<-curr_sero


### Netherlands.
i<-"NLD"
x<-readRDS(paste0("data/derived/",i,"/",i,"_regions.RDS"))
curr_sero<-x$deaths_group
curr_sero$seroprevalence<-x$seroprev_group$seroprevalence
curr_sero$pop<-x$popN*x$prop_pop$pop_prop
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * curr_sero$pop
curr_sero$ifr_crude<-curr_sero$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$deaths_per_million<-1000000*curr_sero$deaths_at_sero/curr_sero$pop
nld_resr<-curr_sero

### Denmark
i<-"DNK"
x<-readRDS(paste0("data/derived/",i,"/",i,"_regions.RDS"))
curr_sero<-x$deaths_group
curr_sero$seroprevalence<-x$seroprev_group$seroprevalence
curr_sero$pop<-x$popN*x$prop_pop$pop_prop
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * curr_sero$pop
curr_sero$ifr_crude<-curr_sero$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$deaths_per_million<-1000000*curr_sero$deaths_at_sero/curr_sero$pop
dnk_resr<-curr_sero

### Switzerland
## process Geneva estimate from the age specific results
x<-readRDS(paste0("data/derived/CHE/CHE_agebands.RDS"))
curr_sero<-x$seroprev
curr_sero$deaths_at_sero<-x$deaths_group$deaths_denom_at_sero[1]
curr_sero$pop<-x$popN
curr_sero$inf_pop_crude<-curr_sero$seroprev * curr_sero$pop
curr_sero$ifr_crude<-curr_sero$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$deaths_per_million<-1000000*curr_sero$deaths_at_sero/curr_sero$pop
che_resr<-curr_sero


## Sweden - to do.

### USA - Los Angeles, California
i<-"LA_CA"
x<-readRDS(paste0("data/derived/USA/",i,"_regions.RDS"))
curr_sero<-x$deaths_group
curr_sero$seroprevalence<-x$seroprev$seroprevalence
curr_sero$pop<-x$popN
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * curr_sero$pop
curr_sero$ifr_crude<-curr_sero$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$deaths_per_million<-1000000*curr_sero$deaths_at_sero/curr_sero$pop
la_ca_resr<-curr_sero

### USA - Santa Clara, California
i<-"SC_CA"
x<-readRDS(paste0("data/derived/USA/",i,"_regions.RDS"))
curr_sero<-x$deaths_group
curr_sero$seroprevalence<-x$seroprev$seroprevalence
curr_sero$pop<-x$popN
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * curr_sero$pop
curr_sero$ifr_crude<-curr_sero$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$deaths_per_million<-1000000*curr_sero$deaths_at_sero/curr_sero$pop
sc_ca_resr<-curr_sero

### USA - CH_MA
i<-"CH_MA"
x<-readRDS(paste0("data/derived/USA/",i,"_regions.RDS"))
curr_sero<-x$deaths_group
curr_sero$seroprevalence<-x$seroprev$seroprevalence
curr_sero$pop<-x$popN
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * curr_sero$pop
curr_sero$ifr_crude<-curr_sero$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$deaths_per_million<-1000000*curr_sero$deaths_at_sero/curr_sero$pop
ch_ma_resr<-curr_sero

### USA - MD_FL
i<-"MD_FL"
x<-readRDS(paste0("data/derived/USA/",i,"_regions.RDS"))
curr_sero<-x$deaths_group
curr_sero$seroprevalence<-x$seroprev$seroprevalence
curr_sero$pop<-x$popN
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * curr_sero$pop
curr_sero$ifr_crude<-curr_sero$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$deaths_per_million<-1000000*curr_sero$deaths_at_sero/curr_sero$pop
md_fl_resr<-curr_sero

### USA - NYC
i<-"NYC_NY_1"
x<-readRDS(paste0("data/derived/USA/",i,"_regions.RDS"))
curr_sero<-x$deaths_group
curr_sero$seroprevalence<-x$seroprev$seroprevalence
curr_sero$pop<-x$popN
curr_sero$inf_pop_crude<-curr_sero$seroprevalence * curr_sero$pop
curr_sero$ifr_crude<-curr_sero$deaths_at_sero/curr_sero$inf_pop_crude
curr_sero$deaths_per_million<-1000000*curr_sero$deaths_at_sero/curr_sero$pop
nyc_resr<-curr_sero



#######################
# PLOTS
######### seroprevalence by age
sero_sheet<-sero_sheet %>%
  dplyr::mutate(n_tested=as.numeric(n_tested))
inds<-which(is.na(sero_sheet$seroprevalence_unadjusted) & is.na(sero_sheet$seroprevalence_weighted))
sero_sheet$seroprevalence_unadjusted[inds]<-sero_sheet$n_positive[inds]/sero_sheet$n_tested[inds]
inds<-which(is.na(sero_sheet$n_positive))
sero_sheet$n_positive[inds]<-sero_sheet$seroprevalence_unadjusted[inds]*sero_sheet$n_tested[inds]

sero_sheetAge<-filter(sero_sheet,age_breakdown==1)
sero_sheetAge<-sero_sheetAge %>%
  dplyr::group_by(study_id,age_low,age_high) %>%
  dplyr::summarise(n_tested=sum(n_tested),
                   n_positive=sum(n_positive),
                   seroprevalence_unadjusted=mean(seroprevalence_unadjusted),
                   seroprevalence_weighted=mean(seroprevalence_weighted),
  )
inds<-which(sero_sheetAge$study_id=="DNK1")
sero_sheetAge$seroprevalence_unadjusted[inds]<-sero_sheetAge$n_positive[inds]/sero_sheetAge$n_tested[inds]
sero_sheetAge$age_high[which(sero_sheetAge$age_high==999)]<-100
sero_sheetAge$age_mid<-0.5*(sero_sheetAge$age_low + sero_sheetAge$age_high)
studies<-unique(sero_sheetAge$study_id)
studies<-c("ESP1","SWE1","CHE1","DNK1","NLD1","GBR2","IRN1","NYC_NY_1", "WENRO_NY_1","LI_NY_1",
           "REST_NY_1")
col_vec <- c(RColorBrewer::brewer.pal(7, "Set1"),cols25(25))
names_studies <- c("Spain", "Sweden", "Switzerland", "Denmark","Netherlands","United Kingdom","Iran",
                   "New York City","WR county, NY","Long Island, NY","Upstate NY")

if(write2file) {
  tiff(file="figures/sero_age.tiff", width=2200,height=1600,res=300,compression="lzw")
  par(mar=c(5,4,4,10))
  plot(sero_sheetAge$age_mid,100*sero_sheetAge$seroprevalence_unadjusted,xlab="age group (years)",ylab="Seroprevalence (%)",
       xlim=c(0,100),col="white")
  for(i in 1:length(studies)) {
    if(i!=6) {
      j<-which(sero_sheetAge$study_id==studies[i])
      points(sero_sheetAge$age_mid[j],100*sero_sheetAge$seroprevalence_unadjusted[j],pch=21, col.main="black", bg=col_vec[i])
      lines(sero_sheetAge$age_mid[j],100*sero_sheetAge$seroprevalence_unadjusted[j],col=col_vec[i])
    }
  }
  legend(105,28,names_studies[c(1:5,7:length(studies))],pch=rep(21,4),col=rep("black",4), bty='n',
         pt.bg=col_vec[c(1:5,7:length(studies))],xpd=T,ncol=1)
  dev.off()
}


######### seroprevalence by gender
sero_sheetGender<-filter(sero_sheet,gender_breakdown==1)
sero_sheetGender<-sero_sheetGender %>%
  dplyr::group_by(study_id,gender) %>%
  dplyr::summarise(n_tested=sum(n_tested),
                   n_positive=sum(n_positive),
                   seroprevalence_unadjusted=mean(seroprevalence_unadjusted),
                   seroprevalence_weighted=mean(seroprevalence_weighted),
                   unique=n()
  )
inds<-which(sero_sheetGender$unique>1)
sero_sheetGender$seroprevalence_unadjusted[inds]<-sero_sheetGender$n_positive[inds]/sero_sheetGender$n_tested[inds]
sero_sheetGender$gender2<-as.factor(sero_sheetGender$gender)
studies<-unique(sero_sheetGender$study_id)
studies<-c("ESP1","SWE1","CHE1","DNK1","NLD1","GBR2","IRN1","NYC_NY_1", "WENRO_NY_1","LI_NY_1",
           "REST_NY_1")
col_vec <- c(RColorBrewer::brewer.pal(7, "Set1"),cols25(25))
names_studies <- c("Spain", "Sweden", "Switzerland", "Denmark","Netherlands","United Kingdom","Iran",
                   "New York City","WR county, NY","Long Island, NY","Upstate NY")

tiff(file="figures/sero_gender.tiff", width=1000,height=1600,res=300,compression="lzw")
par(mar=c(5,4,4,2))
plot(sero_sheetGender$gender2,100*sero_sheetGender$seroprevalence_unadjusted,xlab="gender",ylab="Seroprevalence (%)",
     col="white")
dev.off()


############# DEATHS BY AGE
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
points(la_ca_res$age_mid,scale*la_ca_res$prop_deaths_per_pop,pch=21, col.main="black", bg=col_vec[7])
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
