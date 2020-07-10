library(socialmixr)
source("R/assertions_v5.R")

#' @title Wrapper Function for \link{socialmixr::contact_matrix}
#' @param country string; Name of the country for contact matrix extraction
#' @param surveyDOI string; Survey DOI as returned by \link{socialmixr::get_survey};
#' @param agebands numeric vector; boundaries for age groups
#' @param strict logical; if a country is not included in the survey, should we extrapolate the contact matrix
#' @param nboots_extrapolate numeric; number of boots for extrapolation
#' @details default behavior is to keep missing age groups and sampling from all age groups
#' @details The full listing of survey options from \link{socialmixr::list_surveys} is:
#'          date                                  title             creator                                    url
#'       2020-06-03        Social contact data for Vietnam         Horby Peter https://doi.org/10.5281/zenodo.1289473
#'       2020-06-03           Social contact data for Peru  Carlos G. Grijalva https://doi.org/10.5281/zenodo.1095664
#'       2020-06-03         Social contact data for Russia     Maria Litvinova https://doi.org/10.5281/zenodo.3415222
#'       2020-06-03      Social contact data for Hong Kong        Kathy  Leung https://doi.org/10.5281/zenodo.1165561
#'       2020-06-03             Social contact data for UK Albert Jan van Hoek https://doi.org/10.5281/zenodo.1409506
#'       2020-06-03            POLYMOD social contact data        Joël Mossong https://doi.org/10.5281/zenodo.1043437
#'       2020-06-05 Social contact data for China mainland     Zhang, Juanjuan https://doi.org/10.5281/zenodo.3366396
#'       2020-06-09       Social contact data for Zimbabwe    Alessia Melegaro https://doi.org/10.5281/zenodo.1127693

get_contact_mat <- function(country, surveyDOI, agebands, strict = FALSE, nboots_extrapolate = 5){
  assert_string(country)
  assert_string(surveyDOI)
  assert_numeric(agebands)
  assert_logical(strict)
  if (!strict) {
    assert_numeric(nboots_extrapolate)
  }


  srvy <- socialmixr::get_survey(surveyDOI)
  if(country %in% socialmixr::survey_countries(srvy)) {

    cntmat <- contact_matrix(srvy,
                             countries = country,
                             age.limits = agebands,
                             symmetric = TRUE)

  } else if (!(country %in% socialmixr::survey_countries(srvy)) & !strict) {

    cntmat <- contact_matrix(srvy,
                             countries = socialmixr::survey_countries(srvy),
                             age.limits = agebands,
                             symmetric = TRUE,
                             n = nboots_extrapolate)

    cntmat <- Reduce("+", lapply(cntmat$matrices, function(x) {x$matrix})) / length(cntmat$matrices)

  } else if (!(country %in% socialmixr::survey_countries(srvy)) &  strict) {
    stop("Queried country does not have a contact matrix for the indicated survey")
  }

  return(cntmat)
}

#' @title Helper Function for Regional Contact Counts
get_rgnal_contacts <- function(rgndemog, contactmat) {
  rgndemog.list <- split(rgndemog, factor(rgndemog$region))
  rgnrho.list <- lapply(rgndemog.list, function(x){
    x$popN %*% as.matrix(contactmat)
  })
  rgnrho <- sapply(rgnrho.list, sum)
  return(rgnrho)
}

