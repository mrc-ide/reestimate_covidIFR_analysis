####################################################################################
## Purpose: Recast or aggregate raw data to appropriate marginal distributions
##
## Notes:
####################################################################################
#......................
# BRA1
#......................
BRArgnsero_prevdf <- readr::read_csv("data/raw/Brazil_seroprevalence_first_survey.csv") %>%
  dplyr::select(-c(dplyr::starts_with("X"))) %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(!is.na(positive)) %>% # CHARLIE, note this drop
  dplyr::group_by(region) %>%
  dplyr::summarise(
    n_tested = sum(tested),
    n_positive = sum(positive)
  ) %>%
  dplyr::mutate(
    seroprevalence_unadjusted = n_positive/n_tested,
    country = "BRA",
    study_id = "BRA1",
    age_low = 0,
    age_high = 999,
    gender = "both",
    seroprevalence_weighted = NA,
    date_start_survey = lubridate::ymd("2020-05-15"),
    date_end_survey = lubridate::ymd("2020-05-22"),
    age_breakdown = 0,
    for_regional_analysis = 1,
    gender_breakdown = 0)
write.csv(BRArgnsero_prevdf, "data/raw/RGN_Agg_Brazil_seroprevalence_first_survey.csv")

#......................
# CHE1
#......................
CHE1sero <- readr::read_csv("data/raw/SEROCoV-POP_AggData.csv") %>%
  magrittr::set_colnames(c("studyweek", "gender", "ageband", "positive", "negative", "indeterminate", "dateminobs",	"datemaxobs")) %>%
  dplyr::mutate(date_start_survey = lubridate::dmy(dateminobs),
                date_end_survey = lubridate::dmy(datemaxobs))

CHE1sero_gender <- CHE1sero %>%
  dplyr::group_by(gender, date_start_survey, date_end_survey) %>%
  dplyr::summarise(
    n_positive = sum(positive),
    n_tested = sum(positive) + sum(negative)
  ) %>%
  dplyr::mutate(region = "Geneva",
                age_low = 0,
                age_high = 999,
                age_breakdown = 0,
                gender_breakdown = 1,
                for_regional_analysis = 0,
                seroprevalence_unadjusted = n_positive/n_tested) %>%
  dplyr::select(c("age_low",	"age_high",	"region",	"gender",	"n_tested",	"n_positive", "seroprevalence_unadjusted",
                  "date_start_survey", "date_end_survey", "age_breakdown",	"gender_breakdown", "for_regional_analysis")) %>%
  dplyr::ungroup()

CHE1sero_rgn <- CHE1sero %>%
  dplyr::group_by(date_start_survey, date_end_survey) %>%
  dplyr::summarise(
    n_positive = sum(positive),
    n_tested = sum(positive) + sum(negative)
  ) %>%
  dplyr::mutate(region = "Geneva",
                gender = "both",
                age_low = 0,
                age_high = 999,
                age_breakdown = 0,
                gender_breakdown = 0,
                for_regional_analysis = 1,
                seroprevalence_unadjusted = n_positive/n_tested
                ) %>%
  dplyr::select(c("age_low",	"age_high",	"region",	"gender",	"n_tested",	"n_positive", "seroprevalence_unadjusted",
                  "date_start_survey", "date_end_survey", "age_breakdown",	"gender_breakdown", "for_regional_analysis")) %>%
  dplyr::ungroup()

CHE1sero_ageband <- CHE1sero %>%
  dplyr::mutate(
    age_low = as.numeric( ifelse(ageband == "65+", 65, stringr::str_extract(ageband, "^([0-9]+\\+?)")) ),
    age_high = as.numeric( ifelse(ageband == "65+", 999, stringr::str_extract(ageband, "([0-9]+$)")) )
    ) %>%
  dplyr::group_by(age_low, age_high, date_start_survey, date_end_survey) %>%
  dplyr::summarise(
    n_positive = sum(positive),
    n_tested = sum(positive) + sum(negative)
  ) %>%
  dplyr::mutate(region = "Geneva",
                gender = "both",
                age_breakdown = 1,
                gender_breakdown = 0,
                for_regional_analysis = 0,
                seroprevalence_unadjusted = n_positive/n_tested) %>%
  dplyr::select(c("age_low",	"age_high",	"region",	"gender",	"n_tested",	"n_positive", "seroprevalence_unadjusted",
                  "date_start_survey", "date_end_survey", "age_breakdown",	"gender_breakdown", "for_regional_analysis")) %>%
  dplyr::ungroup()

# out
CHE1sero_full <- dplyr::bind_rows(CHE1sero_gender, CHE1sero_rgn, CHE1sero_ageband)
readr::write_csv(CHE1sero_full, "data/raw/CHE1_seroprev_agg.csv")










#



