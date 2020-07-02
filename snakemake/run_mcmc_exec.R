#!/usr/bin/env Rscript

library(COVIDCurve)
library(drjacoby)
library(optparse)

option_list=list(
  make_option(c("-i", "--inpath"),
              type = "character",
              default = NULL,
              help = paste("input parameters path"),
              metavar = "character"),
  make_option(c("-o", "--outpath"),
              type = "character",
              default = NULL,
              help = paste("out results path"),
              metavar = "character"),
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)







mod <- readRDS(opt$inpath)
#......................
# make cluster object to parallelize chains
#......................
n_chains <- 10
cl <- parallel::makeCluster(n_chains)
fit <- COVIDCurve::run_IFRmodel_agg(IFRmodel = mod$modelobj[[1]],
                                    reparamIFR = TRUE,
                                    reparamInfxn = TRUE,
                                    reparamKnots = TRUE,
                                    chains = n_chains,
                                    burnin = mod$burnin,
                                    samples = mod$samples,
                                    rungs = mod$rungs,
                                    GTI_pow = mod$GTI_pow,
                                    cluster = cl,
                                    silent = FALSE)
mc_accept_mean <- mean(fit$mcmcout$diagnostics$mc_accept$value)
mc_accept_min <- min(fit$mcmcout$diagnostics$mc_accept$value)
# out
out <- list(fit = fit,
            mc_accept_mean = mc_accept_mean,
            mc_accept_min = mc_accept_min)

saveRDS(out, file = opt$outpath)

return(0)

