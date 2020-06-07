library(drake)
source("R/assertions_v5.R")
source("R/process_data.R")
library(drjacoby)
library(COVIDCurve)
plan <- drake_plan(
  #......................
  # assertions that files exist
  #......................

  #......................
  # data input
  #......................
  ESP.agebands.dat <- process_data(deaths = "https://www.dropbox.com/s/z42zfiyr9tr86qe/deaths_v2.csv?dl=1",
                                   population = "https://www.dropbox.com/s/hv4woy1zdlveg72/population.csv?dl=1",
                                   sero_val = "https://www.dropbox.com/s/nu7ek2t2bwo9gxo/seroassay_validation.csv?dl=1",
                                   seroprev = "https://www.dropbox.com/s/kc3mle86os2e6g1/seroprevalence.csv?dl=1",
                                   cumulative = TRUE,
                                   ECDC = "https://www.dropbox.com/s/a2ds6orlpl5ashs/daily_deaths_ECDC20200518.csv?dl=1",
                                   groupingvar = "ageband",
                                   study_ids = "ESP1",
                                   ecdc_countrycode = "ESP",
                                   filtRegions = NULL,
                                   filtGender = NULL,
                                   filtAgeBand = c("0-10", "10-20", "20-30",
                                                   "30-40", "40-50", "50-60",
                                                   "60-70", "70-80", "80-90", "90-999"))
  raw_data = readxl::read_excel(file_in("raw_data.xlsx")),
  data = raw_data %>%
    mutate(Species = forcats::fct_inorder(Species)),
  hist = create_plot(data),
  fit = lm(Sepal.Width ~ Petal.Width + Species, data),
  report = rmarkdown::render(
    knitr_in("report.Rmd"),
    output_file = file_out("report.html"),
    quiet = TRUE
  )
)
plan
