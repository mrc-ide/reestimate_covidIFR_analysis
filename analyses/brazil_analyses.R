# Loading Required Libraries
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Loading Data Dictionary Interconverting Between Regions and States
data_dictionary <- read.csv("C:/Users/cw1716/Documents/COVID_2019/IFR_Update/reestimate_covidIFR_analysis/data/raw/Brazil_State_Region_Data_Dictionary.csv") %>%
  select(-State_Name)

# Loading Brazillian 2020 Population by Age, Sex and State
population <- read.csv("C:/Users/cw1716/Documents/COVID_2019/IFR_Update/reestimate_covidIFR_analysis/data/raw/Brazil_2020_Population_Data.csv")

# Loading Seroprevalence Data
seroprevalence <- read.csv("C:/Users/cw1716/Documents/COVID_2019/IFR_Update/reestimate_covidIFR_analysis/data/raw/Brazil_1st_Seroprevalence_Survey_Results.csv") %>%
  mutate(Positive = case_when(is.na(Positive) ~ Inferred_Positive,
                              TRUE ~ Positive)) %>%
  select(Region, State_Code, Tested, Positive)

# Loading Deaths Data
deaths_data <- readRDS("C:/Users/cw1716/Documents/COVID_2019/IFR_Update/reestimate_covidIFR_analysis/data/derived/BRA/Brazil_state_age_sex_deaths.rds") %>%
  left_join(data_dictionary, by = c("state" = "State"))

# Test Characteristics and Adjustment Function
sens <- 0.864
spec <- 0.996
rogan_gladen <- function(obs_prev, sens, spec) {
  x <- (obs_prev + spec - 1)/(spec + sens - 1)
  return(x)
}

# 1. Calculating National Level IFR
overall_population <- sum(population$Population)
overall_deaths <- sum(deaths_data$count[deaths_data$date <= "2020-05-21"], na.rm = TRUE)
tested <- sum(seroprevalence$Tested)
positive <- sum(seroprevalence$Positive)
overall_seroprevalence <- positive/tested
overall_IFR <- 100 * overall_deaths/(overall_population * overall_seroprevalence)
adj_overall_seroprevalence <- rogan_gladen(overall_seroprevalence, sens, spec)
adj_overall_IFR <- 100 * overall_deaths/(overall_population * adj_overall_seroprevalence)

# 2. Calculating Regional Level IFR
regional_population <- population %>%
  group_by(Region) %>%
  summarise(population = sum(Population))

regional_deaths <- deaths_data %>%
  filter(date <= "2020-05-21") %>%
  group_by(Region) %>%
  summarise(deaths = sum(count))

regional_seroprevalence <- seroprevalence %>%
  group_by(Region) %>%
  summarise(Tested = sum(Tested), Positive = sum(Positive))

overall_region <- regional_seroprevalence %>%
  left_join(regional_deaths, by = "Region") %>%
  left_join(regional_population, by = "Region") %>%
  mutate(seroprevalence = Positive/Tested) %>%
  mutate(IFR = 100 * deaths/(seroprevalence * population)) %>%
  mutate(adj_seroprevalence = rogan_gladen(seroprevalence, sens, spec)) %>%
  mutate(adj_IFR = 100 * deaths/(adj_seroprevalence * population))

# 3. Calculating IFR by Age
age_population <- population %>%
  group_by(Age_Group) %>%
  summarise(Population = sum(Population)) %>%
  spread(Age_Group, Population) %>%
  mutate(`10-20` = `10-15` + `15-20`,
         `20-30` = `20-25` + `25-30`,
         `30-40` = `30-35` + `35-40`,
         `40-50` = `40-45` + `45-50`,
         `50-60` = `50-55` + `55-60`,
         `60-70` = `60-65` + `65-70`,
         `70-80` = `70-75` + `75-80`,
         `80+` = `80-85` + `85-90` + `90+`) %>%
  select(`0-5`, `5-10`, `10-20`, `20-30`, `30-40`, `40-50`, `50-60`, `60-70`, `70-80`, `80+`) %>%
  gather(Age_Group, Population)

# check age bounds are exactly identical to seroprevalence age bounds used
age_deaths_data <- readRDS("C:/Users/cw1716/Documents/COVID_2019/IFR_Update/reestimate_covidIFR_analysis/data/derived/BRA/Brazil_state_age_sex_deaths.rds") %>%
  filter(date <= "2020-05-21") %>%
  mutate(age_category = case_when(
    age >= 0 & age < 5 ~ "0-5",
    age >= 5 & age < 10 ~ "5-10",
    age >= 10 & age < 20 ~ "10-20",
    age >= 20 & age < 30 ~ "20-30",
    age >= 30 & age < 40 ~ "30-40",
    age >= 40 & age < 50 ~ "40-50",
    age >= 50 & age < 60 ~ "50-60",
    age >= 60 & age < 70 ~ "60-70",
    age >= 70 & age < 80 ~ "70-80",
    age >= 80 ~ "80+")) %>%
  group_by(age_category) %>%
  summarise(count = sum(count))

age_seroprevalence <- c(1.4, 1.2, 1.4, 1.4, 1.5, 1.6, 1.7, 1.0, 1.2, 0.5)
age_seroprevalence <- age_seroprevalence/100

overall_age <- age_population %>%
  left_join(age_deaths_data, by = c("Age_Group" = "age_category"))
overall_age$seroprevalence <- age_seroprevalence
overall_age <- overall_age %>%
  mutate(IFR = 100 * count/(seroprevalence * Population)) %>%
  mutate(Age_Group = factor(Age_Group, levels = c("0-5", "5-10", "10-20", "20-30", "30-40",
                                                  "40-50", "50-60", "60-70", "70-80", "80+"))) %>%
  mutate(adj_seroprevalence = rogan_gladen(seroprevalence, sens, spec)) %>%
  mutate(adj_IFR = 100 * count/(adj_seroprevalence * Population))

ggplot(overall_age, aes(x = Age_Group, y = IFR)) +
  geom_bar(stat = "identity") +
  lims(y = c(0, 30))

ggplot(overall_age, aes(x = Age_Group, y = adj_IFR)) +
  geom_bar(stat = "identity") +
  lims(y = c(0, 120))

