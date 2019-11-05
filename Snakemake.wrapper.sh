#!/bin/bash

# to run snakemake as batch job
# run in the data folder for this project

module load snakemake/5.4.4 || exit 1

rm -rf 00log
mkdir -p 00log

sbcmd="sbatch --cpus-per-task={cluster.cpus-per-task} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"


snakefile=$1
config_yaml=$2
cluster_json=$3

snakemake -s $snakefile \
-pr --local-cores 2 --jobs 1999 \
--configfile $config_yaml \
--cluster-config $cluster_json \
--cluster "$sbcmd"  --latency-wait 120  \
-k --restart-times 0

