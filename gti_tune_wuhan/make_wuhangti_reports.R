dir.create(paste0(here::here(), "/gti_tune_wuhan/tune_wuhan_reports/"),
           recursive = TRUE)

#......................
# get results
#......................
paths <- list.files("gti_tune_wuhan/tune_wuhan/", full.names = TRUE)

#......................
# render
#......................
renderMyDocument <- function(path) {
  runchar <- stringr::str_split_fixed(basename(path), "_rung50", n = 2)[[1]]
  rmarkdown::render("gti_tune_wuhan/wuhan_report_base.Rmd",
                    params = list( path = path ),
                    output_file = paste0(here::here(), "/gti_tune_wuhan/tune_wuhan_reports/", runchar, "_wuhangtitune_report.pdf"))
}

lapply(paths, renderMyDocument)


