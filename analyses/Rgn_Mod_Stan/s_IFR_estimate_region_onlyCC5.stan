/////// Seroprevalence regional model for estimating specificity
data {
  int nr;
  int x_seror[nr];
  int N_seror[nr];
  int N_deathsr[nr];
  int tot_obsd;
  vector[nr] popr;
  vector[nr] prop_pop_reg;
  int N_sens_validat;
  int x_sens_validat;
  int N_spec_validat;
  int x_spec_validat;
}

parameters {
  real<lower=0.0001,upper=sum(popr)> tot_inf;
  vector<lower=0.0001,upper=0.9999>[nr] Rne_raw;
  real<lower=0.0001,upper=0.6000> ifr;
  real<lower=0.0001,upper=0.9999> sensitivity;
  real<lower=0.0001,upper=0.9999> specificity;
}

transformed parameters {
  real exptotd;
  real expdr[nr];
  real prev_sero_obsr[nr];
  real prev_sero_truer[nr];
  real exp_infxns_r[nr];
  vector[nr] Rne;

  // reparameterize
  Rne = 0.0001 + Rne_raw .* popr ./ (tot_inf * prop_pop_reg);

  //// Divide up total infections
  for(r in 1:nr) {
      exp_infxns_r[r] = tot_inf * prop_pop_reg[r] *Rne[r];
  }

  for(i in 1:nr) {
    prev_sero_truer[i] = exp_infxns_r[i]/popr[i];
    prev_sero_obsr[i]=(1-specificity)*(1-prev_sero_truer[i]) + sensitivity*prev_sero_truer[i];
  }

    // Marginal deaths
  for(r in 1:nr) {
      expdr[r] = exp_infxns_r[r]*ifr;
  }

  // total deaths likelihood.
  exptotd=sum(expdr);

}

model {
  real prop_exp_dr[nr];

  x_sens_validat ~ binomial(N_sens_validat,sensitivity);
  x_spec_validat ~ binomial(N_spec_validat,specificity);


  for(i in 1:nr) {
    // regional seroprevalence likelihood
    x_seror[i] ~ binomial(N_seror[i],prev_sero_obsr[i]);
  }

  tot_obsd ~ poisson(exptotd);


  for(i in 1:nr) {
    prop_exp_dr[i] = expdr[i]/exptotd;
    // regional deaths likelihood.
    N_deathsr[i] ~ binomial(tot_obsd, prop_exp_dr[i]);
  }
}
