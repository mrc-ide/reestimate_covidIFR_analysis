#!/usr/bin/env bash

# model fits
sbatch -n 1 --time=5-00:00:00 --mem=24g --wrap="Rscript analyses/ModFits/ModFits_drakeworker_v2.R"
sbatch -n 1 --time=5-00:00:00 --mem=24g --wrap="Rscript analyses/ModFits/seroRev_ModFits_drakeworker_v2.R"
sbatch -n 1 --time=5-00:00:00 --mem=18g --wrap="Rscript analyses/ModFits/CareHomes_drakeworker.R"

# simulation fits
sbatch -n 1 --time=8-00:00:00 --mem=24g --wrap="Rscript analyses/SimWork/SimCurves_drakeworker_v2.R"
sbatch -n 1 --time=8-00:00:00 --mem=24g --wrap="Rscript analyses/SimWork/SeroRev_SimCurves_drakeworker_v2.R"
sbatch -n 1 --time=5-00:00:00 --mem=18g --wrap="Rscript analyses/SimWork/fit_conceptual_figure_delayeffects.R"
sbatch -n 1 --time=5-00:00:00 --mem=18g --wrap="Rscript analyses/SimWork/fit_seroday_concept.R"


