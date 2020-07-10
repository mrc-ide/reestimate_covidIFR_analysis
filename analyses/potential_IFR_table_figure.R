library(tidyverse)

#......................
# read in data
#......................
dat <- readxl::read_excel("~/Desktop/small_IFR_table.xlsx") %>%  # temp location
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(!is.na(country))

moddat <- dat %>%
  dplyr::select(c("country", dplyr::starts_with("model"), "population", "seroprevalence")) %>%
  dplyr::mutate(pop = population/1e6)

rtdat <- dat %>%
  dplyr::select(-c(dplyr::starts_with("model"), "population", "seroprevalence")) %>%
  tidyr::gather(., key = "param", value = "est", 2:ncol(.)) %>%
  dplyr::mutate(param = factor(param,
                               levels = c("crude", "capita", "excess"),
                               labels = c("Crude", "Per Capita", "Excess")
                               ))

jpeg("~/Desktop/temp_figure2.jpg", width = 11, height = 8, units = "in", res = 500)
ggplot() +
  geom_point(data = rtdat, aes(x = country, y = est, fill = param, shape = param), size = 6.5, color = "transparent", alpha = 0.8) +
  geom_pointrange(data = moddat,
                  aes(x = country, ymin = model_lci, y = model_median, ymax = model_uci), size = 0.75, alpha = 0.5) +
  geom_point(data = moddat, aes(x = country, y = model_median, size = pop, color = seroprevalence)) +
  ylim(c(0, 3.5)) +
  scale_shape_manual("Unadj. Estimate", values = c(22, 23, 24, 25)) +
  scale_fill_manual("Unadj. Estimate", values = c(wesanderson::wes_palette("IsleofDogs2", type = "discrete"))) +
  scale_size_continuous("Pop. (Millions)", range = c(7, 9), guide = F) + # can't decide if we should keep this or not
  scale_color_gradientn("Mod. IFR (95% CrI) & \n Seroprevalence (%)",
                        colors = rev(wesanderson::wes_palette("Zissou1", 40, type = "continuous")),
                        limits = c(0, 10)) +
  ylab("Infection Fatality Rate") +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  theme(
    plot.title =  element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(family = "Helvetica", face = "bold", vjust = 0.5, hjust = 0.5, size = 16),
    axis.text.x = element_text(family = "Helvetica", color = "#000000", vjust = 0.5, hjust = 0.5, size = 14),
    axis.text.y = element_text(family = "Helvetica", color = "#000000", face = "bold", vjust = 0.5, hjust = 1, size = 14),
    legend.title = element_text(family = "Helvetica", face = "bold", vjust = 0.5, hjust = 0.5, size = 16),
    legend.text = element_text(family = "Helvetica", vjust = 0.8, hjust = 0.5, size = 14, angle = 0),
    legend.position = "right",
    axis.line.x = element_line(color = "black", size = 1.5),
    axis.line.y = element_line(color = "black", size = 1.5),
    axis.ticks.y = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
graphics.off()
