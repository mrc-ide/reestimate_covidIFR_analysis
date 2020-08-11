#......................
# get results
#......................
paths <- list.files("results/ModFits/", full.names = TRUE)

#......................
# render
#......................
renderMyDocument <- function(path) {
  countrychar <- stringr::str_split(basename(path), "_", simplify = T)[[1]]
  stratachar <- ifelse(grepl("age", basename(path)), "Age-Band", "Region")
  rmarkdown::render("reports/report_base.Rmd",
                    params = list( path = path ),
                    output_file = paste0(here::here(), "/reports/", countrychar, "-", stratachar, "_report.pdf"))
}

lapply(paths, renderMyDocument)
