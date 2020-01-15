#!/bin/bash
#SBATCH --time=0-00:15:00
#SBATCH --account=cc-debug
echo Hello, this is job $SLURM_JOB_ID
date
sleep 600
date
