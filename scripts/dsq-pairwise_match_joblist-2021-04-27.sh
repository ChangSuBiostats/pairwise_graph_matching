#!/bin/bash
#SBATCH --array 1-997
#SBATCH --output dsq-pairwise_match_joblist-%A_%3a-%N.out
#SBATCH --job-name dsq-pairwise_match_joblist
#SBATCH --mem-per-cpu 4g -t 4:00 --partition=scavenge
#SBATCH --requeue

# DO NOT EDIT LINE BELOW
/ysm-gpfs/apps/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/ysm/project/fan_zhou/zf59/cs/brain_connectivity/pairwise_graph_matching/scripts/pairwise_match_joblist.txt --status-dir /gpfs/ysm/project/fan_zhou/zf59/cs/brain_connectivity/pairwise_graph_matching/scripts

