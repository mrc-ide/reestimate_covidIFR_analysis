#! /bin/bash

ROOT=/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/ # root directory for project (non-scratch)
SNAKE=/proj/ideel/meshnick/users/NickB/Projects/reestimate_covidIFR_analysis/snakemake
NODES=200 # max number of cluster nodes
WAIT=30 # lag for system

snakemake \
	--snakefile $SNAKE/run_snake_params.py \
	--configfile $SNAKE/config.yaml \
	--directory $ROOT \
	--printshellcmds \
	--rerun-incomplete \
	--keep-going \
	--latency-wait $WAIT \
	--cluster $SNAKE/launch.py \
	-j $NODES \
#	--dryrun -p
