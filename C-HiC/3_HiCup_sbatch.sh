#!/bin/sh

#SBATCH -J ESCs_58
#SBATCH -n 1
#SBATCH -c 6
#SBATCH -p public-cpu
#SBATCH -t 2-00:00:00

srun ${HOME}/Software/hicup_v0.6.1/hicup --config /srv/beegfs/scratch/users/r/roucogar/CHiC_Robert_pipeline/C58ESCs2/1_Config_file_example.conf
