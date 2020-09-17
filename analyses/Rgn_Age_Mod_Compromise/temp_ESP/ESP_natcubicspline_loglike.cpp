Rcpp::List esp_loglike(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc) {

  // extract misc items
  // items needed for death spline
  int n_knots = misc["n_knots"];
  int rcensor_day = misc["rcensor_day"];
  int days_obsd = misc["days_obsd"];
  int agestratlen = misc["agestratlen"];
  int rgnstratlen = misc["rgnstratlen"];

  // items needed serology
  int n_sero_obs = misc["n_sero_obs"];
  std::vector<int> sero_survey_start = Rcpp::as< std::vector<int> >(misc["sero_survey_start"]);
  std::vector<int> sero_survey_end = Rcpp::as< std::vector<int> >(misc["sero_survey_end"]);
  int max_seroday_obsd = misc["max_seroday_obsd"];

  // extract serology items
  double sens = params["sens"];
  double spec = params["spec"];
  double sero_con_rate = params["sero_con_rate"];

  // death delay params
  double mod = params["mod"];
  double sod = params["sod"];

  // extract population demographics by region and put in right format
  std::vector<double> rawdemogpa = Rcpp::as< std::vector<double> >(misc["demog_pa"]);
  std::vector<std::vector<double>> demog_pa(rgnstratlen, std::vector<double>(agestratlen));
  // region A,A,A,  B,B,B  Ages 1,2,3, 1,2,3;
  int iter = 0;
  for (int r = 0; r < rgnstratlen; r++) {
    for (int a = 0; a < agestratlen; a++) {
      demog_pa[r][a] = rawdemogpa[iter];
      iter++;
    }
  }
  // get marg age and regional for populationN
  std::vector<int> Ademog = Rcpp::as< std::vector<int> >(misc["Ademog"]);
  std::vector<int> Rdemog = Rcpp::as< std::vector<int> >(misc["Rdemog"]);

  // storage items
  std::vector<double> ma(agestratlen);
  std::vector<double> Ane(agestratlen);
  std::vector<double> Rne(rgnstratlen);
  std::vector<double> node_x(n_knots);
  std::vector<double> node_y(n_knots);
  // fill storage with parameters
  node_x[0] = 1;
  node_x[1] = params["x1"] *  params["x4"];
  node_x[2] = params["x2"] *  params["x4"];
  node_x[3] = params["x3"] *  params["x4"];
  node_x[4] = params["x4"];
  node_y[0] = params["y1"] * params["y3"];
  node_y[1] = params["y2"] * params["y3"];
  node_y[2] = params["y3"];
  node_y[3] = params["y4"] * params["y3"];
  node_y[4] = params["y5"] * params["y3"];

  ma[0] = params["ma1"] * params["ma10"];
  ma[1] = params["ma2"] * params["ma10"];
  ma[2] = params["ma3"] * params["ma10"];
  ma[3] = params["ma4"] * params["ma10"];
  ma[4] = params["ma5"] * params["ma10"];
  ma[5] = params["ma6"] * params["ma10"];
  ma[6] = params["ma7"] * params["ma10"];
  ma[7] = params["ma8"] * params["ma10"];
  ma[8] = params["ma9"] * params["ma10"];
  ma[9] = params["ma10"];
  // ma[0] = params["ma1"];
  // ma[1] = params["ma2"];
  // ma[2] = params["ma3"];
  // ma[3] = params["ma4"];
  // ma[4] = params["ma5"];
  // ma[5] = params["ma6"];
  // ma[6] = params["ma7"];
  // ma[7] = params["ma8"];
  // ma[8] = params["ma9"];
  // ma[9] = params["ma10"];

  // 10 age groups
  Ane[0] = params["Ane1"];
  Ane[1] = params["Ane2"];
  Ane[2] = params["Ane3"];
  Ane[3] = params["Ane4"];
  Ane[4] = params["Ane5"];
  Ane[5] = params["Ane6"];
  Ane[6] = params["Ane7"];
  Ane[7] = params["Ane8"];
  Ane[8] = params["Ane9"];
  Ane[9] = params["Ane10"];

  // 17 regions in ESP
  Rne[0] = params["Rne1"];
  Rne[1] = params["Rne2"];
  Rne[2] = params["Rne3"];
  Rne[3] = params["Rne4"];
  Rne[4] = params["Rne5"];
  Rne[5] = params["Rne6"];
  Rne[6] = params["Rne7"];
  Rne[7] = params["Rne8"];
  Rne[8] = params["Rne9"];
  Rne[9] = params["Rne10"];
  Rne[10] = params["Rne11"];
  Rne[11] = params["Rne12"];
  Rne[12] = params["Rne13"];
  Rne[13] = params["Rne14"];
  Rne[14] = params["Rne15"];
  Rne[15] = params["Rne16"];
  Rne[16] = params["Rne17"];


  //.............................
  // Assume a uniform attack rate but then allow Ne to rescale
  //.............................
  // ages
  double Anedom = 0.0;
  for (int a = 0; a < agestratlen; a++) {
    Ane[a] = Ane[a] * Ademog[a];
    Anedom += Ane[a];
  }
  // Ne as a prop vector
  for (int a = 0; a < agestratlen; a++) {
    Ane[a] = Ane[a]/Anedom;
  }

  // regions
  double Rnedom = 0.0;
  for (int r = 0; r < rgnstratlen; r++) {
    Rne[r] = Rne[r] * Rdemog[r];
    Rnedom += Rne[r];
  }
  // Ne as a prop vector
  for (int r = 0; r < rgnstratlen; r++) {
    Rne[r] = Rne[r]/Rnedom;
  }

  //.............................
  // Assume age-specific IFR is constant across regions
  //.............................
  std::vector<double> Rma(rgnstratlen);
  // get IFR weighted by region age-demography
  for (int r = 0; r < rgnstratlen; r++) {
    for (int a = 0; a < agestratlen; a++) {
      Rma[r] += ma[a] * demog_pa[r][a];
    }
  }


  //.............................
  // gamma look up table
  //.............................
  std::vector<double> pgmms(days_obsd + 1);
  for (int i = 0; i < (days_obsd+1); i++) {
    pgmms[i] = R::pgamma(i, 1/pow(sod,2), mod*pow(sod,2), true, false);
  }


  //........................................................
  // Deaths Section
  //........................................................
  //.............................
  // knots catch
  //.............................
  const double OVERFLO_DOUBLE = DBL_MAX/100.0;
  double loglik = -OVERFLO_DOUBLE;
  bool nodex_pass = true;
  for (int i = 1; i < node_x.size(); i++) {
    if (node_x[i] <= node_x[i-1]) {
      nodex_pass = false;
    }
  }
  // if (nodex_pass) {
  //.............................
  // natural cubic spline
  //.............................
  // NB, we want N spline functions (interpolants) and so we have n+1 knots
  // as part of simplifying the linear spline, function we write this denom: x_{i+1} - x_i
  std::vector<double> denom(n_knots-1);
  for (int i = 0; i < denom.size(); i++) {
    denom[i] = node_x[i+1] - node_x[i];
  }

  // NB, our intercepts are the y_coordinates of the n+1 knots -- ends of knots will serve as boundaries
  // initialize our knot 2nd derive/linear slope
  std::vector<double> z(n_knots-1);
  z[0] = 0;
  // now store piecewise slopes of second deriv so we can calculate zi's later and make sure we are smoothing through knot
  // -2 b/c first knot and last knot set to 0
  std::vector<double> m(n_knots-2);
  for (int i = 1; i <= m.size(); i++) {
    m[i-1] = (3/denom[i])*(node_y[i+1] - node_y[i]) - (3/denom[i-1])*(node_y[i] - node_y[i-1]);
  }

  // need g and k for first deriv calc of our knots, z
  std::vector<double> g(n_knots-1);
  std::vector<double> k(n_knots-1);
  g[0] = 1;
  k[0] = 0;
  for (int i = 1; i < (n_knots-1); i++) {
    // now we can get our zs and sp2s
    g[i] = 2*(node_x[i+1] - node_x[i-1]) - (denom[i-1])*(k[i-1]);
    k[i] = denom[i]/g[i];
    z[i] = (m[i-1] - denom[i-1]*z[i-1])/g[i];
  }
  // initialize our three "slopes"
  std::vector<double> sp1(n_knots-1);
  std::vector<double> sp3(n_knots-1);
  std::vector<double> sp2(n_knots);
  sp2[n_knots-1] = 0;
  // finally loop through and get our "slopes" for our interpolants
  for (int i = (n_knots-2); i >= 0; i--) {
    sp2[i] = z[i] - k[i]*sp2[i+1];
    sp1[i] = (node_y[i+1] - node_y[i])/(denom[i]) - (denom[i]*(sp2[i+1] + 2*sp2[i]))/3 ;
    sp3[i] = (sp2[i+1] - sp2[i])/(3*denom[i]);
  }

  // create infection spline
  std::vector<double> infxn_spline(days_obsd);
  int node_j = 0;
  infxn_spline[0] = node_y[0];
  for (int i = 1; i < days_obsd; i++) {
    // update curve and have i+1 to account for fact days are 1-based
    infxn_spline[i] = node_y[node_j] +
      sp1[node_j] * ((i+1) - node_x[node_j]) +
      sp2[node_j] * pow(((i+1) - node_x[node_j]), 2) +
      sp3[node_j] * pow(((i+1) - node_x[node_j]), 3);

    // for all interpolants except (potentially) the last knot
    if (node_j < (node_x.size()-2)) {
      // update node_j
      if ((node_x[0] + i) >= node_x[node_j+1]) {
        node_j++;
      }
    }
  }

  // exponentiate infxn spline out of log space
  for (int i = 0; i < days_obsd; i++) {
    infxn_spline[i] = exp(infxn_spline[i]);
  }

  //.............................
  // popN/size catch
  //.............................
  // check if stratified infections exceed stratified population denominator
  bool popN_pass = true;

  std::vector<double> Acum_infxn_check(agestratlen);
  for (int i = 0; i < days_obsd; i++) {
    for (int a = 0; a < agestratlen; a++) {
      Acum_infxn_check[a] += Ane[a] * infxn_spline[i];
    }
  }

  for (int a = 0; a < agestratlen; a++) {
    if (Acum_infxn_check[a] > Ademog[a]) {
      popN_pass = false;
    }
  }

  std::vector<double> Rcum_infxn_check(rgnstratlen);
  for (int i = 0; i < days_obsd; i++) {
    for (int r = 0; r < rgnstratlen; r++) {
      Rcum_infxn_check[r] += Rne[r] * infxn_spline[i];
    }
  }

  for (int r = 0; r < rgnstratlen; r++) {
    if (Rcum_infxn_check[r] > Rdemog[r]) {
      popN_pass = false;
    }
  }

  // if (popN_pass) {

  // loop through days and TOD integral
  std::vector<double> auc(days_obsd);
  for (int i = 0; i < days_obsd; i++) {
    for (int j = i+1; j < (days_obsd + 1); j++) {
      int delta = j - i - 1;
      auc[j-1] += infxn_spline[i] * (pgmms[delta + 1] - pgmms[delta]);
    }
  }

  //.............................
  // L1 Deaths Shape Expectation -- overall
  //.............................
  // items for expectation
  double L1deathshape_loglik = 0.0;
  double texpd = 0.0;
  std::vector<std::vector<double>> agestrata_expdailyd(days_obsd, std::vector<double>(agestratlen));
  std::vector<std::vector<double>> rgnstrata_expdailyd(days_obsd, std::vector<double>(rgnstratlen));



  std::vector<double> expd(days_obsd);
  // read in observed deaths
  std::vector<int> obsd = Rcpp::as< std::vector<int> >(data["obs_deaths"]);
  // get log-likelihood over all days
  for (int  i = 0; i < days_obsd; i++) {

    // get age-specific deaths
    for (int a = 0; a < agestratlen; a++) {
      // store strata deaths for L2
      agestrata_expdailyd[i][a] = auc[i] * Ane[a] * ma[a];
      // get exp deaths per day by "summing out" (Age)strata
      expd[i] += agestrata_expdailyd[i][a];
    }

    // get region-specific deaths
    for (int r = 0; r < rgnstratlen; r++) {
      // store strata deaths for L3
      rgnstrata_expdailyd[i][r] = auc[i] * Rne[r] * Rma[r];
      // note, not summing out deaths again, but by day these should be equal
    }

    // i+1 to account for 1-based dates
    if ((i+1) < rcensor_day) {
      if (obsd[i] != -1) {
        L1deathshape_loglik += R::dpois(obsd[i], expd[i], true);
      }
    }
    // store total expected deaths for L2/L3 likelihood
    texpd += expd[i];
  }

  //.............................
  // L2.1 Deaths Age proportions Expectation
  //.............................
  // items for expectation
  double L2_age_deathprop_loglik = 0.0;
  std::vector<double> agestrata_expd(agestratlen);
  for (int a = 0; a < agestratlen; a++) {
    for (int  i = 0; i < days_obsd; i++) {
      if ((i+1) < rcensor_day) {
        agestrata_expd[a] += agestrata_expdailyd[i][a];
      }
    }
  }

  // extract observed data
  std::vector<double> paobsd = Rcpp::as< std::vector<double> >(data["age_prop_strata_obs_deaths"]);
  for (int a = 0; a < agestratlen; a++) {
    L2_age_deathprop_loglik += R::dbinom(round(agestrata_expd[a]), round(texpd), paobsd[a], true);
  }

  //.............................
  // L2.2 Deaths Rgn proportions Expectation
  //.............................
  // items for expectation
  double L2_rgn_deathprop_loglik = 0.0;
  std::vector<double> rgnstrata_expd(rgnstratlen);
  for (int r = 0; r < rgnstratlen; r++) {
    for (int  i = 0; i < days_obsd; i++) {
      if ((i+1) < rcensor_day) {
        rgnstrata_expd[r] += rgnstrata_expdailyd[i][r];
      }
    }
  }

  // extract observed data
  std::vector<double> probsd = Rcpp::as< std::vector<double> >(data["rgn_prop_strata_obs_deaths"]);
  for (int r = 0; r < rgnstratlen; r++) {
    L2_rgn_deathprop_loglik += R::dbinom(round(rgnstrata_expd[r]), round(texpd), probsd[r], true);
  }

  //........................................................
  // Serology Section
  //........................................................
  // get cumulative hazard for each day up to the latest serology observation date
  // i.e. cumulative hazard of seroconversion on given day look up table
  std::vector<double> cum_serocon_hazard(max_seroday_obsd);
  for (int d = 0; d < max_seroday_obsd; d++) {
    cum_serocon_hazard[d] = 1-exp((-(d+1)/sero_con_rate));
  }

  // seroconversion by strata look up table
  std::vector<std::vector<double>> sero_con_num_full(max_seroday_obsd, std::vector<double>(rgnstratlen));
  // loop through and split infection curve by strata and by number of seroconversion study dates
  for (int r = 0; r < rgnstratlen; r++) {
    for (int i = 0; i < max_seroday_obsd; i++) {
      // go to the "end" of the day
      for (int j = i+1; j < (max_seroday_obsd + 1); j++) {
        int time_elapsed = j - i - 1;
        sero_con_num_full[j-1][r] += infxn_spline[i] * Rne[r] * cum_serocon_hazard[time_elapsed];
        // sero_con_num_full[j-1][r] -= infxn_spline[i] * Rne[r] * cum_serorev_hazard[time_elapsed];
      }
    }
  }

  // get average over serostudy data
  std::vector<std::vector<double>> sero_con_num(n_sero_obs, std::vector<double>(rgnstratlen));
  for (int i = 0; i < n_sero_obs; i++) {
    for (int j = 0; j < rgnstratlen; j++) {
      for (int k = sero_survey_start[i]; k <= sero_survey_end[i]; k++) {
        // days are 1 based
        sero_con_num[i][j] += sero_con_num_full[k-1][j];
      }
      // now get average
      sero_con_num[i][j] =  sero_con_num[i][j]/(sero_survey_end[i] - sero_survey_start[i] + 1);
    }
  }

  //.............................
  // L3.1 AGE Serology Proportions Expectation
  //.............................
  double L3_agesero_loglik = 0.0;
  // unpack serology observed data
  std::vector<double> age_datpos_raw = Rcpp::as< std::vector<double> >(data["age_obs_serologypos"]);
  std::vector<double> age_datn_raw = Rcpp::as< std::vector<double> >(data["age_obs_serologyn"]);
  // recast datpos
  std::vector<std::vector<double>> age_datpos(n_sero_obs, std::vector<double>(agestratlen));
  std::vector<std::vector<double>> age_datn(n_sero_obs, std::vector<double>(agestratlen));
  int ageseroiter = 0;
  for (int i = 0; i < n_sero_obs; i++) {
    for (int a = 0; a < agestratlen; a++) {
      age_datpos[i][a] = age_datpos_raw[ageseroiter];
      age_datn[i][a] = age_datn_raw[ageseroiter];
      ageseroiter++;
    }
  }
  // loop through sero likelihood
  for (int i = 0; i < n_sero_obs; i++) {
    for (int a = 0; a < agestratlen; a++) {
      if (age_datpos[i][a] != -1 | age_datn[i][a] != -1 ) {
        // Gelman Estimator for numerical stability
        double obs_prev = sens*(sero_con_num[i][a]/Ademog[a]) + (1-spec)*(1 - (sero_con_num[i][a]/Ademog[a]));
        L3_agesero_loglik += R::dbinom(age_datpos[i][a], age_datn[i][a], obs_prev, true);
      }
    }
  }


  //.............................
  // L3.1 RGN Serology Proportions Expectation
  //.............................
  double L3_rgnsero_loglik = 0.0;
  // unpack serology observed data
  std::vector<double> rgn_datpos_raw = Rcpp::as< std::vector<double> >(data["rgn_obs_serologypos"]);
  std::vector<double> rgn_datn_raw = Rcpp::as< std::vector<double> >(data["rgn_obs_serologyn"]);
  // recast datpos
  std::vector<std::vector<double>> rgn_datpos(n_sero_obs, std::vector<double>(rgnstratlen));
  std::vector<std::vector<double>> rgn_datn(n_sero_obs, std::vector<double>(rgnstratlen));
  int rgnseroiter = 0;
  for (int i = 0; i < n_sero_obs; i++) {
    for (int r = 0; r < rgnstratlen; r++) {
      rgn_datpos[i][r] = rgn_datpos_raw[rgnseroiter];
      rgn_datn[i][r] = rgn_datn_raw[rgnseroiter];
      rgnseroiter++;
    }
  }
  // loop through sero likelihood
  for (int i = 0; i < n_sero_obs; i++) {
    for (int r = 0; r < rgnstratlen; r++) {
      if (rgn_datpos[i][r] != -1 | rgn_datn[i][r] != -1 ) {
        // Gelman Estimator for numerical stability
        double obs_prev = sens*(sero_con_num[i][r]/Rdemog[r]) + (1-spec)*(1 - (sero_con_num[i][r]/Rdemog[r]));
        L3_rgnsero_loglik += R::dbinom(rgn_datpos[i][r], rgn_datn[i][r], obs_prev, true);
      }
    }
  }


  // bring together
  loglik = L1deathshape_loglik + L2_age_deathprop_loglik + L2_rgn_deathprop_loglik + L3_agesero_loglik + L3_rgnsero_loglik;

  // catch underflow
  if (!std::isfinite(loglik)) {
    loglik = -OVERFLO_DOUBLE;
  }

  // end cumulative vs. popN check
  // }
  // end node_x check
  // }
  // return loglike
  // return Rcpp::wrap(loglik);
  Rcpp::List ret = Rcpp::List::create(L1deathshape_loglik, L2_age_deathprop_loglik, L2_rgn_deathprop_loglik, L3_agesero_loglik, L3_rgnsero_loglik);
  return ret;
}
