#!/bin/bash
#
#SBATCH --account=freedman_mlf3 #note: this may need to be updated depending on the user
#SBATCH --job-name=chips
#SBATCH --output=out.chips.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=medium
#SBATCH --mem=30G
#SBATCH -t 2-0
#SBATCH --propagate=STACK


ulimit -s unlimited #added 5/12/21 to fix a problem with segmentation fault error with sambamba due to exceeding stack limit

srun --account=freedman_mlf3 snakemake -s chips/chips.snakefile -j 100 -pr --rerun-incomplete --latency-wait 60 --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {threads} --account=freedman_mlf3"



