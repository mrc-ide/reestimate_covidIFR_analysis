## .................................................................................
## Purpose: Compare meta-analysis IFR estimates
##
## Notes:
## .................................................................................
library(tidyverse)
source("R/covidcurve_helper_functions.R")
source("R/my_themes.R")
#......................
# our data
#......................
brzetal <- readr::read_tsv("tables/final_tables/overall_best_est_for_age_IFRs.tsv") %>%
  dplyr::mutate(age = purrr::map_dbl(ageband, get_mid_age)) %>%
  tidyr::pivot_longer(data = ., cols = -c("ageband", "age"), names_to = "modlvl", values_to = "est") %>%
  dplyr::mutate(model = ifelse(grepl("reg", modlvl), "Brazeau_etal_NoSerorev", "Brazeau_etal_Serorev"),
                Q2.5 = as.numeric(stringr::str_extract(est, "[0-9]+\\.?[0-9]+?(?=\\,)")),
                Q97.5 = as.numeric(stringr::str_extract(est, "[0-9]+\\.?[0-9]+?(?=\\))")),
                Q50 = as.numeric(stringr::str_extract(est, "[0-9]+\\.?[0-9]+?(?=\\s)")),
                # catch non decimals
                Q2.5 = ifelse(is.na(Q2.5), as.numeric(stringr::str_extract(est, "[0-9]+(?=\\,)")), Q2.5),
                Q97.5 = ifelse(is.na(Q97.5), as.numeric(stringr::str_extract(est, "[0-9]+(?=\\)")), Q97.5),
                Q50 = ifelse(is.na(Q50), as.numeric(stringr::str_extract(est, "[0-9]+(?=\\s)")), Q50),
                ) %>%
  dplyr::select(c("age", "model", "Q2.5", "Q50", "Q97.5"))

# other studies
ifrstudies <- readr::read_tsv("meta_comparison/updated_exterior_ifr_table.tsv") %>%
  dplyr::mutate(age = purrr::map_dbl(age, function(x){
    low <- as.numeric(stringr::str_split_fixed(x, "-", n = 2)[,1])
    high <- as.numeric(stringr::str_split_fixed(x, "-", n = 2)[,2]) + 1
    high <- ifelse(high == 1000, 100, high)
    return((low+high)/2)
  })) %>%
  dplyr::bind_rows(., brzetal)


#......................
# plot out
#......................
ifrstudies %>%
  ggplot() +
  geom_pointrange(aes(x = age, y = Q50, ymin = Q2.5, ymax = Q97.5,
                      color = model), alpha = 0.8) +
  scale_color_viridis_d("Study") +
  ylab("IFR (95% Interval") + xlab("Age") +
  #scale_y_log10() +
  xyaxis_plot_theme +
  theme(legend.position = "bottom")
