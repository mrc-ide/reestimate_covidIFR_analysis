#......................
# get results assuming "seroinfinity"
#......................
paths <- list.files("results/Modfits_noserorev/", full.names = TRUE)

#......................
# render
#......................
renderMyDocument <- function(path) {
  countrychar <- stringr::str_split(basename(path), "_", simplify = T)[[1]]
  rmarkdown::render("reports/Modfit_report_base.Rmd",
                    params = list( path = path ),
                    output_file = paste0(here::here(), "/reports/", countrychar, "-", "_NoSeroRev_report.pdf"))
}

lapply(paths, renderMyDocument)



#......................
# get results assuming seroreversion
#......................
paths <- list.files("results/Modfits_serorev/", full.names = TRUE)

#......................
# render
#......................
renderMyDocument <- function(path) {
  countrychar <- stringr::str_split(basename(path), "_", simplify = T)[[1]]
  rmarkdown::render("reports/Modfit_report_base.Rmd",
                    params = list( path = path ),
                    output_file = paste0(here::here(), "/reports/", countrychar, "-", "_SeroRev_report.pdf"))
}

lapply(paths, renderMyDocument)

#......................
# get results for excluding care home deaths
#......................
paths <- list.files("results/Modfits_carehomes/", full.names = TRUE)

#......................
# render
#......................
renderMyDocument <- function(path) {
  countrychar <- stringr::str_split(basename(path), "_", simplify = T)[[1]]
  rmarkdown::render("reports/Modfit_report_base.Rmd",
                    params = list( path = path ),
                    output_file = paste0(here::here(), "/reports/", countrychar, "-", "_CareHomes_report.pdf"))
}

lapply(paths, renderMyDocument)


