dir.create(paste0(here::here(), "/reports/sims/"), recursive = TRUE)
#......................
# get results assuming "seroinfinity"
#......................
paths <- list.files("results/SimCurves_noserorev/", full.names = TRUE)

#......................
# render
#......................
renderMyDocument <- function(path) {
  countrychar <- stringr::str_split(basename(path), "_", simplify = T)[[1]]
  stratachar <- ifelse(grepl("age", basename(path)), "Age-Band", "Region")
  rmarkdown::render("reports/Sim_report_base.Rmd",
                    params = list( path = path ),
                    output_file = paste0(here::here(), "/reports/sims/", countrychar, "-", stratachar, "_NoSeroRev_report.pdf"))
}

lapply(paths, renderMyDocument)



#......................
# get results assuming seroreversion
#......................
paths <- list.files("results/SimCurves_serorev/", full.names = TRUE))

#......................
# render
#......................
renderMyDocument <- function(path) {
  countrychar <- stringr::str_split(basename(path), "_", simplify = T)[[1]]
  stratachar <- ifelse(grepl("age", basename(path)), "Age-Band", "Region")
  rmarkdown::render("reports/Sim_report_base.Rmd",
                    params = list( path = path ),
                    output_file = paste0(here::here(), "/reports/sims/", countrychar, "-", stratachar, "_SeroRev_report.pdf"))
}

lapply(paths, renderMyDocument)



