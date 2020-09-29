#### Regional data
library(tidyverse)
regdf<-read.csv("data/raw/jrc-covid-19-all-days-by-regions.csv")
sero_prevdf <- readr::read_tsv("data/raw/seroprevalence_final_raw.tsv") %>%
  dplyr::select(-c("ref", "notes")) %>%
  dplyr::mutate(date_start_survey = lubridate::ymd(date_start_survey), # NB, we just convert this to a lubridate format and later within the process data function, dates are converted to international format
                date_end_survey = lubridate::ymd(date_end_survey))

regdf<-regdf %>%
  dplyr::mutate(Date = lubridate::ymd(Date)) %>%
  dplyr::filter(CountryName %in% c("Spain","United Kingdom","Italy"))
readr::write_csv(regdf,"data/raw/jrc-covid-19-all-days-by-regions.csv")

#plot(regdf$Date,regdf$CumulativeDeceased)

#### Spain
# sero timings
serotimes<-sero_prevdf %>%
  dplyr::filter(country=="Spain") %>%
  dplyr::select(date_start_survey,date_end_survey) %>%
  dplyr::filter(!duplicated(date_start_survey,date_end_survey)) %>%
  rowwise %>%
  dplyr::mutate(seromid=mean.Date(c(date_start_survey,date_end_survey))) %>%
  dplyr::pull(seromid)

regdf %>%
  dplyr::filter(CountryName=="Spain") %>%
  ggplot( aes(x=Date, y=CumulativeDeceased, group=Region, color=Region)) +
  geom_line() + theme_bw() + scale_y_log10() +
  geom_vline(xintercept=serotimes,linetype="dashed")

## Italy
serotimes<-sero_prevdf %>%
  dplyr::filter(country=="Italy") %>%
  dplyr::select(date_start_survey,date_end_survey) %>%
  dplyr::filter(!duplicated(date_start_survey,date_end_survey)) %>%
  rowwise %>%
  dplyr::mutate(seromid=mean.Date(c(date_start_survey,date_end_survey))) %>%
  dplyr::pull(seromid)

regdf %>%
  dplyr::filter(CountryName=="Italy") %>%
  ggplot( aes(x=Date, y=CumulativeDeceased, group=Region, color=Region)) +
  geom_line() + theme_bw() + scale_y_log10() +
  geom_vline(xintercept=serotimes,linetype="dashed")

## UK - missing deaths?
serotimes<-sero_prevdf %>%
  dplyr::filter(study_id=="GBR3") %>%
  dplyr::select(date_start_survey,date_end_survey) %>%
  dplyr::filter(!duplicated(date_start_survey,date_end_survey)) %>%
  rowwise %>%
  dplyr::mutate(seromid=mean.Date(c(date_start_survey,date_end_survey))) %>%
  dplyr::pull(seromid)

england<-c("London","North West","Midlands","North East and Yorkshire","South East","South West",
           "East of England")
regdf %>%
  dplyr::filter(CountryName=="United Kingdom" & Region %in% england)
%>%
  ggplot( aes(x=Date, y=CumulativeDeceased, group=Region, color=Region)) +
  geom_line() + theme_bw() + #scale_y_log10() +
  geom_vline(xintercept=serotimes,linetype="dashed")
