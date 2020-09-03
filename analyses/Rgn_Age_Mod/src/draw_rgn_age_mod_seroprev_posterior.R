#............................................................
#
#...........................................................
source("R/assertions_v5.R")
RgnAgeMod_draw_posterior_sero_curves <- function(rgn_age_IFRmodObj, whichrung = "rung1", dwnsmpl, by_chain = TRUE) {
  assert_pos_int(dwnsmpl)
  assert_string(whichrung)
  assert_logical(by_chain)
  #......................
  # fitler to sampling and by rung
  #......................
  mcmcout.nodes <- rgn_age_IFRmodObj$mcmcout$output %>%
    dplyr::filter(stage == "sampling" & rung == whichrung)

  #......................
  # sample by CI limit and make infxn curves
  #......................
  mcmcout.nodes <- mcmcout.nodes %>%
    dplyr::mutate(logposterior = loglikelihood + logprior)
  # Log-Sum-Exp trick
  convert_post_probs <- function(logpost) {
    exp(logpost - (log(sum(exp(logpost - max(logpost)))) + max(logpost)))
  }
  probs <- convert_post_probs(mcmcout.nodes$logposterior)
  # downsample
  dwnsmpl_rows <- sample(1:nrow(mcmcout.nodes), size = dwnsmpl,
                         prob = probs)
  dwnsmpl_rows <- sort(dwnsmpl_rows)
  mcmcout.nodes <- mcmcout.nodes[dwnsmpl_rows, ]

  #............................................................
  # cpp manipulation
  #...........................................................
  draw_posterior_seroprev <- readLines("analyses/Rgn_Age_Mod/temp_ESP/ESP_natcubicspline_loglike.cpp")
  # remove comments which cause issues when coercing to string in this format
  commlines <- grep("//", draw_posterior_seroprev)
  draw_posterior_seroprev <- draw_posterior_seroprev[! 1:length(draw_posterior_seroprev) %in% commlines]
  draw_posterior_seroprev <- gsub("\\\n", "", draw_posterior_seroprev)
  draw_posterior_seroprev <- capture.output(cat(draw_posterior_seroprev))
  draw_posterior_seroprev <- sub("SEXP", "Rcpp::List", draw_posterior_seroprev)
  draw_posterior_seroprev <- stringr::str_replace(draw_posterior_seroprev,"if \\(nodex_pass\\) \\{", "")
  draw_posterior_seroprev <- stringr::str_replace(draw_posterior_seroprev, "  if \\(popN_pass\\) \\{", "")
  draw_posterior_seroprev <- stringr::str_split_fixed(draw_posterior_seroprev, "std::vector<double> cum_serocon_hazard\\(max_seroday_obsd\\);", n = 2)[,1]

  # rewriting the sero_con_num_full vector here to be all days observed, not just serology days
  draw_posterior_seroprev <- paste(draw_posterior_seroprev,
                                   "std::vector<double> cum_serocon_hazard(days_obsd);
                           for (int d = 0; d < days_obsd; d++) {
                             cum_serocon_hazard[d] = 1-exp((-(d+1)/sero_con_rate));
                           }
                           std::vector<std::vector<double>> sero_con_num_full(days_obsd, std::vector<double>(rgnstratlen));
                            for (int r = 0; r < rgnstratlen; r++) {
                              for (int i = 0; i < days_obsd; i++) {
                                for (int j = i+1; j < (days_obsd + 1); j++) {
                                  int time_elapsed = j - i - 1;
                                  sero_con_num_full[j-1][r] += infxn_spline[i] * Rne[r] * cum_serocon_hazard[time_elapsed];
                                }
                              }
                            }
                           std::vector<std::vector<double>> crude_seroprev(days_obsd, std::vector<double>(rgnstratlen));
                           std::vector<std::vector<double>> RG_seroprev(days_obsd, std::vector<double>(rgnstratlen));
                            for (int i = 0; i < days_obsd; i++) {
                              for (int j = 0; j < rgnstratlen; j++) {
                                crude_seroprev[i][j] = (sero_con_num_full[i][j]/demog_margrgn[j]);
                                RG_seroprev[i][j] = sens*crude_seroprev[i][j] + (1-spec)*(1-crude_seroprev[i][j]);
                              }
                            }",
                                   "Rcpp::List ret = Rcpp::List::create(sero_con_num_full, crude_seroprev,  RG_seroprev); return ret;}",
                                   collapse = "")
  Rcpp::cppFunction(draw_posterior_seroprev)

  #......................
  # inputs needed for cpp function
  #......................
  # misc list
  misc_list = rgn_age_IFRmodObj$misc_list
  # data in
  datin <- rgn_age_IFRmodObj$datinput

  #......................
  # split, run, recombine
  #......................
  cpp_function_wrapper <- function(params, datin, misc) {
    paramsin <- params[!names(params) %in% c("chain", "rung", "iteration", "stage", "logprior", "loglikelihood", "logposterior")]
    paramsin <- unlist(paramsin)

    seroprev_lists <- esp_loglike(params = paramsin,
                                  param_i = 1,
                                  data = datin,
                                  misc = misc_list)

    # extract relevant bits
    sero_counts <- seroprev_lists[[1]] %>%
      do.call("rbind.data.frame", .) %>%
      magrittr::set_colnames(paste0("serocounts_", paste0("R", 1:sum(grepl("Rne", names(paramsin)))))) %>%
      dplyr::mutate(ObsDay = 1:length(datin$obs_deaths)) %>%
      dplyr::select(c("ObsDay", dplyr::everything()))

    crude_seroprev <- seroprev_lists[[2]] %>%
      do.call("rbind.data.frame", .) %>%
      magrittr::set_colnames(paste0("crude_pd_seroprev_",paste0("R", 1:sum(grepl("Rne", names(paramsin)))))) %>%
      dplyr::mutate(ObsDay = 1:length(datin$obs_deaths)) %>%
      dplyr::select(c("ObsDay", dplyr::everything()))

    RG_seroprev <- seroprev_lists[[3]] %>%
      do.call("rbind.data.frame", .) %>%
      magrittr::set_colnames(paste0("RG_pd_seroprev_", paste0("R", 1:sum(grepl("Rne", names(paramsin)))))) %>%
      dplyr::mutate(ObsDay = 1:length(datin$obs_deaths)) %>%
      dplyr::select(c("ObsDay", dplyr::everything()))

    # out
    ret <- dplyr::left_join(sero_counts, crude_seroprev, by = "ObsDay") %>%
      dplyr::left_join(., RG_seroprev, by = "ObsDay")
    return(ret)

  }

  mcmcout.node.rows <- split(mcmcout.nodes, 1:nrow(mcmcout.nodes))
  mcmcout.nodes$seroprev <- purrr::map(mcmcout.node.rows, cpp_function_wrapper,
                                       datin = datin, misc = misc_list)

  #......................
  # tidy
  #......................
  # keep params around for convenience
  if (by_chain) {
    dat <- mcmcout.nodes %>%
      dplyr::select(c("chain", "seroprev")) %>%
      dplyr::group_by(chain) %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      dplyr::ungroup(chain) %>%
      tidyr::unnest(cols = "seroprev")


  } else {
    # keep params around for convenience
    dat <- mcmcout.nodes %>%
      dplyr::select("seroprev") %>%
      dplyr::mutate(sim = 1:dplyr::n()) %>%
      tidyr::unnest(cols = "seroprev")
  }

  #......................
  # out
  #......................
  return(dat)
}
