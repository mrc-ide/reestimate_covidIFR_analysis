data {
  int nr;
  int na;    // number of age groups for IFR calculation
  int na_s;        // number of sero survey age groups (not always the same as deaths)
  int d_i[na];      // index for how deaths by age map to serology age groups.
  int x_seror[nr];
  int N_seror[nr];
  int x_seroa[na_s];
  int N_seroa[na_s];
  int N_deathsr[nr];
  int N_deathsa[na];
  int tot_obsd;
  vector[nr] popr;
  vector[nr] prop_pop_reg;
  vector[na] popa;
  vector[na] prop_pop_age;
  vector[na_s] popas;// population in age groups for the serosurvey.
  vector[na_s] prop_pop_as;
  int pop_reg_age[nr,na];
  real prop_pop_reg_age[nr,na];
  real prop_pop_reg_age_sero[nr,na_s];
  int N_sens_validat;
  int x_sens_validat;
  int N_spec_validat;
  int x_spec_validat;
}

parameters {
  real<lower=0.0001,upper=sum(popr)> tot_inf;
  vector<lower=0.0001,upper=0.9999>[nr] Rne_raw;
  vector<lower=0.0001,upper=0.9999>[na_s] Ane_raw;

  real<lower=0.0001,upper=0.40> ifr[na];
  real<lower=0.0001,upper=0.9999> sensitivity;
  real<lower=0.0001,upper=0.9999> specificity;
}

transformed parameters {
  real exptotd;
  real expdr[nr];
  real expda[na];
  real prev_sero_obsr[nr];
  real prev_sero_obsa[na_s];
  real prev_sero_truer[nr];
  real prev_sero_truea[na_s];
  real exp_infxns_r[nr];
  real exp_infxns_a[na];
  real exp_infxns_as[na_s];
  real exp_infxns[nr,na];
  real exp_infxns_ras[nr,na_s];
  vector[nr] Rne;
  vector[na_s] Ane;

  // rescale params to distribute infections by region and age, acc to upper limit in each group
  Rne = 0.0001 + Rne_raw .* popr ./ (tot_inf * prop_pop_reg);
  Ane = 0.0001 + Ane_raw .* popas ./ (tot_inf * prop_pop_as);


  //// Divide up total infections into regions and serological age groups.
  for(r in 1:nr) {
    for(a_s in 1:na_s) {
      exp_infxns_ras[r,a_s] = tot_inf * prop_pop_reg_age_sero[r,a_s] *Rne[r]*Ane[a_s];
    }
  }

  //// Further divide infections amongst deaths age groups.
  for(r in 1:nr) {
    for(a in 1:na) {
      exp_infxns[r,a] = exp_infxns_ras[r,d_i[a]] * prop_pop_reg_age[r,a]/prop_pop_reg_age_sero[r,d_i[a]];
    }
  }

  ///// Add up infections by age and region
  for(r in 1:nr) exp_infxns_r[r] = 0.0;
  for(a in 1:na) exp_infxns_a[a] =0.0;
  for(r in 1:nr) {
    for(a in 1:na) {
      exp_infxns_r[r] +=exp_infxns[r,a];
      exp_infxns_a[a] +=exp_infxns[r,a];
    }
  }

  /// Seroprevalence likelihood
  for(i in 1:nr) {
    prev_sero_truer[i] = exp_infxns_r[i]/popr[i];
    prev_sero_obsr[i]=(1-specificity)*(1-prev_sero_truer[i]) + sensitivity*prev_sero_truer[i];
  }
  /// Add up infections by age strata in sero data
  for(a_s in 1:na_s) exp_infxns_as[a_s]=0.0;

  for(r in 1:nr) {
    for(a_s in 1:na_s) {
      exp_infxns_as[a_s] += exp_infxns_ras[r,a_s];
    }
  }

  for(i in 1:na_s) {
    prev_sero_truea[i] = exp_infxns_as[i]/popas[i];
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

  // total deaths
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

  for(i in 1:na_s) {
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
