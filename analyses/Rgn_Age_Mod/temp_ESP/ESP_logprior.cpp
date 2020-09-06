SEXP logprior(Rcpp::NumericVector params, int param_i, Rcpp::List misc) {
  double ma1 = params["ma1"];
  double ma2 = params["ma2"];
  double ma3 = params["ma3"];
  double ma4 = params["ma4"];
  double ma5 = params["ma5"];
  double ma6 = params["ma6"];
  double ma7 = params["ma7"];
  double ma8 = params["ma8"];
  double ma9 = params["ma9"];
  double ma10 = params["ma10"];
  double x1 = params["x1"];
  double x2 = params["x2"];
  double x3 = params["x3"];
  double x4 = params["x4"];
  double y1 = params["y1"];
  double y2 = params["y2"];
  double y3 = params["y3"];
  double y4 = params["y4"];
  double y5 = params["y5"];
  double sens = params["sens"];
  double spec = params["spec"];
  double sero_con_rate = params["sero_con_rate"];
  double Rne1 = params["Rne1"];
  double Rne2 = params["Rne2"];
  double Rne3 = params["Rne3"];
  double Rne4 = params["Rne4"];
  double Rne5 = params["Rne5"];
  double Rne6 = params["Rne6"];
  double Rne7 = params["Rne7"];
  double Rne8 = params["Rne8"];
  double Rne9 = params["Rne9"];
  double Rne10 = params["Rne10"];
  double Rne11 = params["Rne11"];
  double Rne12 = params["Rne12"];
  double Rne13 = params["Rne13"];
  double Rne14 = params["Rne14"];
  double Rne15 = params["Rne15"];
  double Rne16 = params["Rne16"];
  double Rne17 = params["Rne17"];
  double mod = params["mod"];
  double sod = params["sod"];
  double ret = R::dunif(ma1,0,1,true) +
    R::dunif(ma2,0,1,true) +
    R::dunif(ma3,0,1,true) +
    R::dunif(ma4,0,1,true) +
    R::dunif(ma5,0,1,true) +
    R::dunif(ma6,0,1,true) +
    R::dunif(ma7,0,1,true) +
    R::dunif(ma8,0,1,true) +
    R::dunif(ma9,0,1,true) +
    R::dunif(ma10,0,1,true) +
    R::dunif(x1,0,0.33,true) +
    R::dunif(x2,0.33,0.66,true) +
    R::dunif(x3,0.66,0.99,true) +
    R::dunif(x4,175,200,true) +
    R::dunif(y1,0,1,true) +
    R::dunif(y2,0,1,true) +
    R::dunif(y3,0,1,true) +
    R::dunif(y4,0,1,true) +
    R::dunif(y5,0,10,true) +
    R::dnorm(sero_con_rate,18,1, true) +
    R::dbeta(sens,123.5,30.5, true) +
    R::dbeta(spec,156.5,0.5, true) +
    R::dnorm(Rne1,1,0.05,true) +
    R::dnorm(Rne2,1,0.05,true) +
    R::dnorm(Rne3,1,0.05,true) +
    R::dnorm(Rne4,1,0.05,true) +
    R::dnorm(Rne5,1,0.05,true) +
    R::dnorm(Rne6,1,0.05,true) +
    R::dnorm(Rne7,1,0.05,true) +
    R::dnorm(Rne8,1,0.05,true) +
    R::dnorm(Rne9,1,0.05,true) +
    R::dnorm(Rne10,1,0.05,true) +
    R::dnorm(Rne11,1,0.05,true) +
    R::dnorm(Rne12,1,0.05,true) +
    R::dnorm(Rne13,1,0.05,true) +
    R::dnorm(Rne14,1,0.05,true) +
    R::dnorm(Rne15,1,0.05,true) +
    R::dnorm(Rne16,1,0.05,true) +
    R::dnorm(Rne17,1,0.05,true) +
    R::dnorm(mod,19.5,1,true) +
    R::dbeta(sod,79,21,true) +
    3*log(x4) +
    4*log(y3);
  if (!std::isfinite(ret)) {
    const double OVERFLO_DOUBLE = DBL_MAX/100.0;
    ret = -OVERFLO_DOUBLE;
  }
  return Rcpp::wrap(ret);
}



















//
