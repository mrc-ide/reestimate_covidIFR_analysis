/////// Seroprevalence regional model for estimating specificity
data {
  int nr;
  int na;
  int x_seror[nr];
  int N_seror[nr];
  int x_seroa[na];
  int N_seroa[na];
  int N_deathsr[nr];
  int N_deathsa[na];
  int tot_obsd;
  vector[nr] popr;
  vector[nr] prop_pop_reg;
  vector[na] popa;
  vector[na] prop_pop_age;
  int pop_reg_age[nr,na];
  real prop_pop_reg_age[nr,na];
  int N_sens_validat;
  int x_sens_validat;
  int N_spec_validat;
  int x_spec_validat;
}

parameters {
  real<lower=0.0001,upper=sum(popr)> tot_inf;
  vector<lower=0.0001,upper=0.9999>[nr] Rne_raw;
  vector<lower=0.0001,upper=0.9999>[na] Ane_raw;
  real<lower=0.0001,upper=0.9999> ifr[na];
  real<lower=0.0001,upper=0.9999> sensitivity;
  real<lower=0.0001,upper=0.9999> specificity;
}

transformed parameters {
  real exptotd;
  real expdr[nr];
  real expda[na];
  real prev_sero_obsr[nr];
  real prev_sero_obsa[na];
  real<lower=0.0001,upper=0.9999> prev_sero_truer[nr];
  real prev_sero_truea[na];
  real exp_infxns_r[nr];
  real exp_infxns_a[na];
  real exp_infxns[nr,na];
  vector[nr] Rne;
  vector[na] Ane;
  
  Rne = 0.0001 + Rne_raw .* popr ./ (tot_inf * prop_pop_reg);
  Ane = 0.0001 + Ane_raw .* popa ./ (tot_inf * prop_pop_age);
  
  //// Divide up total infections
  for(r in 1:nr) {
    for(a in 1:na) {
      exp_infxns[r,a] = tot_inf * prop_pop_reg_age[r,a] *Rne[r]*Ane[a];
    } 
  }
  
  //// Marginal true seroprev
  // add up infections
  for(r in 1:nr) exp_infxns_r[r] = 0.0;
  for(a in 1:na) exp_infxns_a[a] =0.0;
  for(r in 1:nr) {
    for(a in 1:na) {
      exp_infxns_r[r] +=exp_infxns[r,a];
      exp_infxns_a[a] +=exp_infxns[r,a];
    }
  }

  for(i in 1:nr) {
    prev_sero_truer[i] = exp_infxns_r[i]/popr[i];
    prev_sero_obsr[i]=(1-specificity)*(1-prev_sero_truer[i]) + sensitivity*prev_sero_truer[i];
  }
  for(i in 1:na) {
    prev_sero_truea[i] = exp_infxns_a[i]/popa[i];
    prev_sero_obsa[i]=(1-specificity)*(1-prev_sero_truea[i]) + sensitivity*prev_sero_truea[i];
  }
    // Marginal deaths
    for(r in 1:nr) expdr[r]=0.0;
  for(r in 1:nr) {
    for(a in 1:na) {
      expdr[r] += exp_infxns[r,a]*ifr[a];
    }
  }

  for(i in 1:na) expda[i] = 0.0;
  for(a in 1:na) {
    for(r in 1:nr) {
      expda[a] += exp_infxns[r,a]*ifr[a];
    }
  }

  // total deaths likelihood.
  exptotd=sum(expdr);

}

model {
  real prop_exp_dr[nr];
  real prop_exp_da[na];
  
  x_sens_validat ~ binomial(N_sens_validat,sensitivity);
  x_spec_validat ~ binomial(N_spec_validat,specificity);
  
  
  for(i in 1:nr) {
    // regional seroprevalence likelihood
    x_seror[i] ~ binomial(N_seror[i],prev_sero_obsr[i]);
  }
  
  for(i in 1:na) {
    // age seroprevalence likelihood
    x_seroa[i] ~ binomial(N_seroa[i],prev_sero_obsa[i]);
  }

  tot_obsd ~ poisson(exptotd);
  
  
  for(i in 1:nr) {
    prop_exp_dr[i] = expdr[i]/exptotd;
    // regional deaths likelihood.
    N_deathsr[i] ~ binomial(tot_obsd, prop_exp_dr[i]);  
  }
  for(i in 1:na) {
    prop_exp_da[i] = expda[i]/exptotd;
    // age deaths likelihood.
    N_deathsa[i] ~ binomial(tot_obsd, prop_exp_da[i]);  
  }
}
