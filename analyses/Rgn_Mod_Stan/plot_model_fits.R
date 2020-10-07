
library(rstan)

#/// MODEL FITS REGIONAL MODEL - PLOTS

studies<-c("ESP1-2","ITA1","GBR3","DNK1","NYS1","BRA1")

for(i in 1:length(studies)) {
  curr_study_id<-studies[i]
  fit<-readRDS(paste0("C:/Users/Lucy/Dropbox (SPH Imperial College)/IFR update/rgn_mod_results/final_fits/fit_",curr_study_id,"_reg_age_full.rds"))
  dat<-readRDS(paste0("analyses/Rgn_Mod_Stan/input_data/",curr_study_id,"input_dat.rds"))
  params<-rstan::extract(fit)
  prev_sero_reg<-100*apply(params$prev_sero_truer,2,mean)
  prev_sero_reg_lci<-100*apply(params$prev_sero_truer,2,quantile,0.025)
  prev_sero_reg_uci<-100*apply(params$prev_sero_truer,2,quantile,0.975)
  prev_sero_reg_obs<-100*apply(params$prev_sero_obsr,2,mean)
  prev_sero_reg_obs_lci<-100*apply(params$prev_sero_obsr,2,quantile,0.025)
  prev_sero_reg_obs_uci<-100*apply(params$prev_sero_obsr,2,quantile,0.975)
  if(!is.null(params$prev_sero_truea)) {
    prev_sero_age<-apply(params$prev_sero_truea,2,mean)
    prev_sero_age_obs<-apply(params$prev_sero_obsa ,2,mean)
    prev_sero_age_obs_uci<-apply(params$prev_sero_obsa,2,quantile,0.975)
  }
  expdr<-100000*apply(params$expdr ,2,mean)/dat$popr
  expdr_lci<-100000*apply(params$expdr ,2,quantile,0.025)/dat$popr
  expdr_uci<-100000*apply(params$expdr ,2,quantile,0.975)/dat$popr


  plot2file<-T
  if(plot2file) { # overwrite
    tiff(paste0("analyses/Rgn_Mod_Stan/results/" ,curr_study_id, "_regional_fit.tiff"),width=3000,height=950,compression="lzw",res=300)
    layout(matrix(c(1:4), nrow = 1, ncol = 4, byrow = TRUE), widths=c(2,2,1.5,2))

    max_x<-0.5+max(prev_sero_reg_uci,100*dat$x_seror/dat$N_seror)
    plot(100*dat$x_seror/dat$N_seror,100000*dat$N_deathsr/dat$popr,pch=19,
         ylab="deaths per 100,000",xlab="Seroprevalence (%)",
         ylim=c(0,max(expdr_uci)),
         xlim=xlim)
    plotCI(prev_sero_reg,expdr,
           li=prev_sero_reg_lci,
           ui=prev_sero_reg_uci,
           col="cornflowerblue",err="x",add=T,sfrac=0,gap=0.01)
    plotCI(prev_sero_reg,expdr,
           li=expdr_lci,
           ui=expdr_uci,
           col="cornflowerblue", err="y",add=T,sfrac=0,gap=0.01)
    legend(1,max(expdr_uci)*1.3,c("data","fitted (test-adjusted)"),pch=c(19,1),col=c("black","cornflowerblue"),xpd=T,bty='n')
    legend(-max_x/13,max(expdr_uci),"A",bty="n")

    ########## Seroprev (fitted) vs mortality
    plot(100*dat$x_seror/dat$N_seror,100000*dat$N_deathsr/dat$popr,pch=19,
         ylab="deaths per 100,000",xlab="Seroprevalence (%)",
         ylim=c(0,max(expdr_uci,100*dat$x_seror/dat$N_seror)),
         xlim=c(0,0.5+max(prev_sero_reg_obs_uci,100*dat$x_seror/dat$N_seror)))
    plotCI(prev_sero_reg_obs,expdr,
           li=prev_sero_reg_obs_lci,
           ui=prev_sero_reg_obs_uci,
           add=T,col="firebrick",err="x",cex=0.8,sfrac=0,gap=0.01)
    plotCI(prev_sero_reg_obs,expdr,
           li=expdr_lci,
           ui=expdr_uci,
           col="firebrick", err="y",add=T,sfrac=0,gap=0.01)
    legend(1,max(expdr_uci)*1.3,c("data","fitted (observed)"),pch=c(19,1),col=c("black","firebrick"),xpd=T,bty='n')
    legend(-max_x/13,max(expdr_uci),"B",bty="n")

    ############ Sens and spec
    d<-density(params$specificity)
    plot(d$x,d$y/sum(d$y),xlim=c(0.95,1.005),type="l",xlab="Specificity",ylab="density",main="",lwd=0.9)
    lines(rep(dat$x_spec_validat/dat$N_spec_validat,2),c(0,500),lty=2)
    legend(0.943,max(d$y/sum(d$y)),"C",bty="n")


    ########## Seroprev age
    if(curr_study_id!="DNK1") {
      y1<-dat$x_seroa/dat$N_seroa
      sens<-mean(params$sensitivity)
      spec<-mean(params$specificity)
      y3<-prev_sero_age*sens + (1-prev_sero_age)*(1-spec)

      max_x<-length(y1)+0.5
      age_labels<-dat$sero_agebands
      age_labels<-gsub("-999","+",age_labels)
      plot(1:length(y1)-0.1, y1, pch=19, col="black",add=T,ylim=c(0,max(prev_sero_age_obs_uci+0.05)),xaxt='n',ylab="Seroprevalence",
           xlab="",xlim=c(0.5,max_x))
      plotCI(1:length(y1)+0.1,prev_sero_age_obs,
             li=apply(params$prev_sero_obsa,2,quantile,0.025),
             ui=prev_sero_age_obs_uci,
             col="cornflowerblue",add=T,pch=19,sfrac=0)
      axis(1,at=1:length(y1),labels = age_labels,las=2,)
      mtext("Age",side=1,line=3.5,cex=0.65)
      legend(1,max(prev_sero_age_obs_uci+0.05)*1.3,c("data","fitted (observed)"),pch=c(19,19),col=c("black","cornflowerblue"),xpd=T,bty='n')
      legend((max_x-0.5)/20+0.01,max(prev_sero_age_obs_uci+0.05),"D",bty="n")
    }


    dev.off()

  }

}


####################
# schematic of relationship of seroprevalence and mortality with different test sensitivity and specificity.
####################
if(plot2file) {
  jpeg("analyses/Rgn_Mod_Stan/results/regional_schematic.jpg",width=1600,height=1300,res=300)
  true_sero<-1:30/100
  sens<-0.8
  spec<-0.98
  obs_sero<-true_sero*sens + (1-true_sero)*(1-spec)
  deaths<-true_sero*0.007*100000

  par(mar=c(5,4,5,2))
  plot(true_sero*100,deaths,type="l",xlab="Seroprevalence % (observed or true)", ylab="deaths per 100,000",lwd=2)
  lines(obs_sero*100,deaths,col="blue")
  sens<-1
  obs_sero<-true_sero*sens + (1-true_sero)*(1-spec)
  lines(obs_sero*100,deaths,col="red")
  sens<-0.8
  spec<-1
  obs_sero<-true_sero*sens + (1-true_sero)*(1-spec)
  lines(obs_sero*100,deaths,col="orange")
  legend(0,300, c("true relationship","observed, sens=80%, spec=98%", "observed, sens=100%, spec=98%",
                  "observed, sens=80%, spec=100%"),
         col=c("black","blue","red","orange"),lty=1,bty='n',xpd=T,cex=0.85)
  dev.off()
}


### PLOT SENSITIVITY
# d<-density(params$sensitivity)
# plot(d$x,d$y/sum(d$y),xlim=c(0.7,1.005),type="l",xlab="parameter",ylab="density",main="")
# lines(rep(dat$x_sens_validat/dat$N_sens_validat,2),c(0,500),lty=2)
#
