#! /usr/bin/env python3
"""
Master script for running slurm workers
"""

from __future__ import print_function

import os
import sys
import yaml
import re
from collections import defaultdict

def load_run_metadata(f):
	""" Get params from tab-separated file."""
	params = list()
	with open(f) as rfile:
		for line in rfile:
			if line.startswith("#"):
				continue
			line = line.strip().split("\t")
			params.append(line[0])
	return params


## read run manifest and neccessary supports
parampath = load_run_metadata(config["manifest"])
param_in_dir = config["paramindir"]
param_out_dir = config["outdir"]

final_target = []
for i in parampath:
	final_target.append( os.path.join(param_out_dir, "{}_results.RDS".format(i)) )


rule all:
	input: final_target

rule params_out:
	input: param_in_dir + "{params}.RDS"
	output: param_out_dir + "{params}_results.RDS",
	log: param_out_dir + "{params}_log.Rout",
	shell:
		r"""
		Rscript --max-ppsize=500000 --vanilla \
			snakemake/run_mcmc_exec.R \
			--inpath {input} \
			--outpath {output} \
			>& {log}
		"""
