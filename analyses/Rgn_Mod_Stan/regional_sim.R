
setwd("C:/Users/Lucy/Documents/GitHub/reestimate_covidIFR_analysis")

library(rstan)
library(shinystan)
library(tidyr)
library(plotrix)
library(reshape)

### Load STAN model
model_reg_only_full<-stan_model("analyses/Rgn_Mod_Stan/s_IFR_estimate_region_onlyCC5.stan")

#################################
# Test performance & example validation study
#################################
spec_true<-0.985
sens_true<-0.8
N_sens_validat<-100   
N_spec_validat<-100
x_sens_validat<-80      # test sensitivity is 80%
x_spec_validat<-100   #fix to be wrong. rbinom(1,N_spec_validat,spec_true)
sens<-x_sens_validat/N_sens_validat  ## observed sens and spec in validation study
spec<-x_spec_validat/N_spec_validat


#################################
# Simulate regional seroprevalence and deaths
#################################
nr<-10 # number of regions

## population
tot_pop<-5800000   ## total population across regions
probs_reg<-exp(runif(nr,log(0.01),log(1))) # generate uneven probability of being in each region - varying population size
pop_reg<-rmultinom(1,tot_pop,probs_reg)  ## assign population to regions


### seroprevalence
## generate a range of prevalences up to 10%. Most regions experience smaller epidemics.
prev_sero_true<-exp(seq(log(0.001),log(0.1),length.out=nr))   #sample on log scale.
prev_sero_obs_prob<-prev_sero_true*sens_true + (1-prev_sero_true)*(1-spec_true) #generate observed prevalence

## deaths
ifr<-0.007
mean_deaths<-pop_reg*prev_sero_true*ifr
N_deaths_reg<- rpois(nr,mean_deaths)
N_sero_tot<-60000 # total sample size for serosurvey.
N_sero_reg<-rmultinom(1,N_sero_tot,probs_reg)  ## serosample in proportion to region size
x_sero_reg<-apply(cbind(N_sero_reg,prev_sero_obs_prob),1, function(x) rbinom(1,x[1],x[2]))

#plot(x_sero_reg/N_sero_reg,N_deaths_reg/pop_reg)

sim_data<-data.frame(x_sero_reg=x_sero_reg,N_sero_reg=N_sero_reg,
                     N_deaths_reg=N_deaths_reg,pop_reg=pop_reg)
saveRDS(sim_data, file="analyses/Rgn_Mod_Stan/sim_reg_dat_spec.rds")


nIter<-10000

fit_reg_age_full <- sampling(model_reg_only_full,list(nr=nr,
                                                      x_seror=x_sero_reg,
                                                      N_seror=as.vector(N_sero_reg),
                                                      N_deathsr=N_deaths_reg,
                                                      tot_obsd = sum(N_deaths_reg),
                                                      popr=as.vector(pop_reg),
                                                      prop_pop_reg = as.vector(pop_reg/sum(pop_reg)),
                                                      x_sens_validat=x_sens_validat,
                                                      N_sens_validat=N_sens_validat,
                                                      x_spec_validat=x_spec_validat,
                                                      N_spec_validat=N_spec_validat),
                             iter=nIter, chains=4,control = list(adapt_delta = 0.95, max_treedepth = 10))

params<-rstan::extract(fit_reg_age_full)
prev_sero_reg<-100*apply(params$prev_sero_truer,2,mean)
prev_sero_reg_lci<-100*apply(params$prev_sero_truer,2,quantile,0.025)
prev_sero_reg_uci<-100*apply(params$prev_sero_truer,2,quantile,0.975)
prev_sero_reg_obs<-100*apply(params$prev_sero_obsr,2,mean)
prev_sero_reg_obs_lci<-100*apply(params$prev_sero_obsr,2,quantile,0.025)
prev_sero_reg_obs_uci<-100*apply(params$prev_sero_obsr,2,quantile,0.975)

expdr<-100000*apply(params$expdr ,2,mean)/pop_reg
expdr_lci<-100000*apply(params$expdr ,2,quantile,0.025)/pop_reg
expdr_uci<-100000*apply(params$expdr ,2,quantile,0.975)/pop_reg




####################################
##### PLOT TOGETHER WITH SCHEMATIC

plot2file<-T
if(plot2file) { # overwrite
  #jpeg(paste0("analyses/Rgn_Mod_Stan/results/" ,curr_study_id, "_regional_fit.jpg"),width=3000,height=950,res=300)
  #layout(matrix(c(1:4), nrow = 1, ncol = 4, byrow = TRUE), widths=c(2,2,1.5,2))
  jpeg("analyses/Rgn_Mod_Stan/results/regional_fit_sim.jpg",width=2500,height=950,res=300)

  
  par(mfrow=c(1,3))  
  true_sero<-0:30/100
  sens_true<-0.8
  spec_true<-0.985
  obs_sero<-true_sero*sens_true + (1-true_sero)*(1-spec_true)
  deaths<-true_sero*0.007*100000
  
  par(mar=c(5,4,5,2))
  plot(true_sero*100,deaths,type="l",xlab="seroprevalence (%)", ylab="deaths per 100,000",lwd=2)
  lines(obs_sero*100,deaths,col="blue")
  sens_true<-1
  obs_sero<-true_sero*sens_true + (1-true_sero)*(1-spec_true)
  lines(obs_sero*100,deaths,col="red")
  sens_true<-0.8
  spec_true<-1
  obs_sero<-true_sero*sens_true + (1-true_sero)*(1-spec_true)
  lines(obs_sero*100,deaths,col="orange")
  legend(-3,300, c("true","observed, imperfect sens & spec", "observed, perfect sens, imperfect spec",
                  "observed, perfect spec, imperfect sens"),
         col=c("black","blue","red","orange"),lty=1,bty='n',xpd=T)
  arrows(1.5, -45, y1=-10, length = 0.1, xpd = T)
  mtext("(1-specificity)",line=2.2,side=1,at=1.5,cex=0.6)
  mtext("(A)", side=2,at=270,line=1.5,las=1,font=2)
  
  sens_true<-0.8
  spec_true<-0.985
  max_x<-0.5+max(prev_sero_reg_uci,100*x_sero_reg/N_sero_reg)
  plot(100*x_sero_reg/N_sero_reg,100000*N_deaths_reg/pop_reg,pch=19,
       ylab="deaths per 100,000",xlab="seroprevalence (%)",
       ylim=c(0,max(expdr_uci)),
       xlim=c(0,max_x))
  plotCI(prev_sero_reg,expdr,
         li=prev_sero_reg_lci,
         ui=prev_sero_reg_uci,
         col="cornflowerblue",err="x",add=T,sfrac=0,gap=0,pch=19)
  plotCI(prev_sero_reg,expdr,
         li=expdr_lci,
         ui=expdr_uci,
         col="cornflowerblue", err="y",add=T,sfrac=0,gap=0)
  legend(1,max(expdr_uci)*1.3,c("simulated data","fitted (test-adjusted)"),pch=19,col=c("black","cornflowerblue"),xpd=T,bty='n')
  mtext("(B)", side=2,at=98,line=1.5,las=1,font=2)
  
  
  ############ Inferred posterior spec
  d<-density(params$specificity)
  plot(d$x*100,d$y/sum(d$y),xlim=c(97,100.5),type="l",xlab="specificity",ylab="density",main="",lwd=0.9)
  lines(rep(100*x_spec_validat/N_spec_validat,2),c(0,500),lty=2)
  lines(rep(spec_true*100,2),c(0,500),lty=2,col="blue")
  mtext("(C)", side=2,at=0.0083,line=1.5,las=1,font=2)
  legend(97,0.009, c("initial spec estimate","true spec", 
                  "posterior spec"),
         col=c("black","blue","black"),lty=c(2,2,1),bty='n',xpd=T)
  
  dev.off()  
}