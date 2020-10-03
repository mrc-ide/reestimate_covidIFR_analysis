####################################################################################
## Purpose: Basic themes for IFR project
##
## Notes:
####################################################################################

# white background but strong x-y axes
xyaxis_plot_theme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 14),
                          axis.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, size = 12),
                          axis.text = element_text(family = "Helvetica", hjust = 0.5, size = 11),
                          legend.position = "right",
                          legend.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5, vjust = 0.5, size = 12),
                          legend.background = element_blank(),
                          legend.text = element_text(family = "Helvetica", hjust = 0, vjust = 0.5, size = 10),
                          legend.key = element_blank(),
                          panel.background = element_rect(fill = "transparent"),
                          plot.background = element_rect(fill = "transparent"),
                          panel.grid = element_blank(),
                          panel.border = element_blank(),
                          axis.line = element_line(color = "#000000", size = 1))
