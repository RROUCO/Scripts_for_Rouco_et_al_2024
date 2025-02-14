#!/bin/sh

module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load SAMtools/1.9
module load R/3.5.1

salloc -n1 -c1 --partition=shared-cpu srun sh /srv/beegfs/scratch/users/r/roucogar/CHiC_Robert_pipeline/CFL58DPrep3E115/Generate_profiles_FLDPE115.sh