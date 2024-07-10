#!/bin/bash
#SBATCH --partition=norm
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=8g
#SBATCH --time=96:00:00
#SBATCH --no-requeue

#export PATH=/mnt/nasapps/development/python/3.7.1/bin:$PATH
#export PATH=/mnt/ccrsf-ifx/Software/tools/python/3.7.1_oel8/bin:$PATH
module load snakemake/5
snakemake --jobname 's.{jobid}.{rulename}'  --latency-wait 600 -k --stats stats/snakemake_$(date +"%Y%m%d%H%M%S").stats --rerun-incomplete --restart-times 0 -j 300 --printshellcmds    --cluster-config config/cluster.yaml  --cluster "sbatch --gres {cluster.gres} --cpus-per-task {cluster.threads} -p {cluster.partition} -t {cluster.time} --mem {cluster.mem} --job-name {cluster.name} --output {cluster.output} --error {cluster.error}"  >&  logs/snakemake_$(date +"%Y%m%d%H%M%S").log

