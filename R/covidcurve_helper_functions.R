source("R/assertions_v5.R")
library(tidyverse)
library(stringr)

#' @title Make Simple Data Dictionary Key for IFR age-bands, regions, etc. to simple
#' @param strata string vector; Names of stata to simplify

make_ma_dict_key <- function(strata_names) {
  assert_string(strata_names)
  tibble::tibble(strata_name = strata_names,
                 param_name  = paste0("ma", 1:length(strata_name)))
}


#' @title Make IFR Uniform Distributed Reparameterized Param Df
#' @param num_mas positive interger; Number of IFR strata to infer

make_ma_reparamdf <- function(num_mas = 10) {
  assert_pos_int(num_mas)
  tibble::tibble(name = paste0("ma", 1:num_mas),
                 min  = rep(0, num_mas),
                 init = rep(0.5, num_mas),
                 max = rep(1, num_mas),
                 dsc1 = rep(0, num_mas),
                 dsc2 = rep(1, num_mas))
}


#' @title Make Spline Y Position Reparameterized Param Df
#' @param num_ys positive interger; Number of Spline Y positions to infer
#' @param max_yvec vector; Describes the ymax position. Rest will be uniform 0,1 for reparameterized positions

make_spliney_reparamdf <- function(max_yvec = list("name" = "y3", min = 0, init = 9, max = 14, dsc1 = 0, dsc2 = 14),
                              num_ys = 5) {
  assert_pos_int(num_ys)
  assert_in(names(max_yvec), c("name", "min", "init", "max", "dsc1", "dsc2"))
  assert_in(c("name", "min", "init", "max", "dsc1", "dsc2"), names(max_yvec))
  assert_string(max_yvec[["name"]])
  assert_numeric(max_yvec[["min"]])
  assert_numeric(max_yvec[["init"]])
  assert_numeric(max_yvec[["max"]])
  assert_numeric(max_yvec[["dsc1"]])
  assert_numeric(max_yvec[["dsc2"]])

  out <- tibble::tibble(name = paste0("y", 1:num_ys),
                        min  = rep(0, num_ys),
                        init = rep(0.5, num_ys),
                        max = rep(1, num_ys),
                        dsc1 = rep(0, num_ys),
                        dsc2 = rep(1, num_ys))
  out %>%
    dplyr::filter(name != max_yvec["name"]) %>%
    dplyr::bind_rows(., max_yvec) %>%
    dplyr::arrange(name)
}

#' @title Make Spline X (Knots) Position Reparameterized Param Df
#' @param num_xs positive interger; Number of Spline X positions to infer
#' @param max_xvec vector; Describes the ymax position. Rest will be uniform 0,1 for reparameterized positions

make_splinex_reparamdf <- function(max_xvec = list("name" = "x4", min = 180, init = 190, max = 200, dsc1 = 180, dsc2 = 200),
                                   num_xs = 4) {
  assert_pos_int(num_xs)
  assert_in(names(max_xvec), c("name", "min", "init", "max", "dsc1", "dsc2"))
  assert_in(c("name", "min", "init", "max", "dsc1", "dsc2"), names(max_xvec))
  assert_string(max_xvec[["name"]])
  assert_numeric(max_xvec[["min"]])
  assert_numeric(max_xvec[["init"]])
  assert_numeric(max_xvec[["max"]])
  assert_numeric(max_xvec[["dsc1"]])
  assert_numeric(max_xvec[["dsc2"]])

  out <- tibble::tibble(name = paste0("x", 1:num_xs),
                        min  = rep(0, num_xs),
                        init = rep(0.5, num_xs),
                        max = rep(1, num_xs),
                        dsc1 = rep(0, num_xs),
                        dsc2 = rep(1, num_xs))
  out %>%
    dplyr::filter(name != max_xvec["name"]) %>%
    dplyr::bind_rows(., max_xvec) %>%
    dplyr::arrange(name)
}

