
library(rstan)

#/// MODEL FITS REGIONAL MODEL - PLOTS
studies<-c("ESP1-2","ITA1","GBR3","DNK1","NYS1","BRA1")
curr_study_id<-studies[4]
fit<-readRDS(paste0("C:/Users/Lucy/Dropbox (SPH Imperial College)/IFR update/rgn_mod_results/final_fits/fit_reg_age_full_",curr_study_id,".rds"))
dat<-readRDS(paste0("analyses/Rgn_Mod_Stan/input_data/",curr_study_id,"input_dat.rds"))
params<-rstan::extract(fit)
prev_sero_reg<-100*apply(params$prev_sero_truer,2,mean)
prev_sero_reg_lci<-100*apply(params$prev_sero_truer,2,quantile,0.025)
prev_sero_reg_uci<-100*apply(params$prev_sero_truer,2,quantile,0.975)
prev_sero_reg_obs<-100*apply(params$prev_sero_obsr,2,mean)
prev_sero_reg_obs_lci<-100*apply(params$prev_sero_obsr,2,quantile,0.025)
prev_sero_reg_obs_uci<-100*apply(params$prev_sero_obsr,2,quantile,0.975)
prev_sero_age<-apply(params$prev_sero_truea,2,mean)
prev_sero_age_obs<-apply(params$prev_sero_obsa ,2,mean)
expdr<-100000*apply(params$expdr ,2,mean)/dat$popr
expdr_lci<-100000*apply(params$expdr ,2,quantile,0.025)/dat$popr
expdr_uci<-100000*apply(params$expdr ,2,quantile,0.975)/dat$popr


par(mfrow=c(1,1))
plotCI(prev_sero_reg,expdr,
       li=prev_sero_reg_lci,
       ui=prev_sero_reg_uci,
       col="lightblue",err="x",ylab="deaths per 100,000",xlab="Seroprevalence (%)",
       ylim=c(0,max(expdr_uci)),
       xlim=c(0,max(prev_sero_reg_uci)),sfrac=0)
plotCI(prev_sero_reg,expdr,
       li=expdr_lci,
       ui=expdr_uci,
       col="lightblue",pch=19, err="y",add=T,sfrac=0)
points(100*dat$x_seror/dat$N_seror,100000*dat$N_deathsr/dat$popr,pch=19)

########## Seroprev (fitted) vs mortality
par(mfrow=c(1,1))
plotCI(prev_sero_reg_obs,expdr,
       li=prev_sero_reg_obs_lci,
       ui=prev_sero_reg_obs_uci,
       col="blue",err="x",ylab="deaths per 100,000",xlab="Seroprevalence (%)",
       ylim=c(0,max(expdr_uci)),
       xlim=c(0,max(c(prev_sero_reg_obs_uci,100*x_sero_reg/N_sero_reg))),sfrac=0)
plotCI(prev_sero_reg_obs,expdr,
       li=expdr_lci,
       ui=expdr_uci,
       col="blue",pch=19, err="y",add=T,sfrac=0)
plotCI(prev_sero_reg,expdr,
       li=prev_sero_reg_lci,
       ui=prev_sero_reg_uci,
       col="orange",err="x",ylab="deaths per 100,000",xlab="Seroprevalence (%)",
       ylim=c(0,max(100000*N_deaths_reg/pop_reg)),
       xlim=c(0,max(prev_sero_reg_uci)),pch=19,sfrac=0,add=T)
points(100*dat$x_seror/dat$N_seror,100000*dat$N_deathsr/dat$popr,pch=19)

########## Seroprev age
y1<-dat$x_seroa/dat$N_seroa
sens<-mean(params$sensitivity)
spec<-mean(params$specificity)
y3<-prev_sero_age*sens + (1-prev_sero_age)*(1-spec)

age_labels<-as.character(unique(dat_age$seroprev_group$ageband))
age_labels<-gsub("-999","+",age_labels)
plot(1:length(y1), y1, pch=19, col="black",add=T,ylim=c(0,max(y1+0.05)),xaxt='n',ylab="Observed seroprevalence",
     xlab="")
plotCI(1:length(y1),prev_sero_age_obs,
       li=apply(params$prev_sero_obsa,2,quantile,0.025),
       ui=apply(params$prev_sero_obsa,2,quantile,0.975),
       col="blue",add=T,pch=19,sfrac=0)
axis(1,at=1:length(y1),labels = age_labels,las=2)

############ Sens and spec
d<-density(params$specificity)
plot(d$x,d$y/sum(d$y),xlim=c(0.95,1.005),type="l",xlab="parameter",ylab="density",main="")
lines(rep(dat$x_spec_validat/dat$N_spec_validat,2),c(0,500),lty=2)

### PLOT SENSITIVITY
d<-density(params$sensitivity)
plot(d$x,d$y/sum(d$y),xlim=c(0.7,1.005),type="l",xlab="parameter",ylab="density",main="")
lines(rep(dat$x_sens_validat/dat$N_sens_validat,2),c(0,500),lty=2)

