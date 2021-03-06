#!/bin/bash
#PBS -N jz_cov_omp
##PBS -l nodes=node2067.grid4cern.if.pw.edu.pl:ppn=40
##PBS -l nodes=1:ppn=40
#PBS -l walltime=4:00:00
#PBS -j oe
## [NAME_OF_THE_JOB].o [NAME_OF_THE_JOB].e


## ------ QUEUE SYSTEM ------
## to submit use:
##     qsub job.st.dwarf
## to check status use:
##     qstat 
## to delete job use:
##     qdel jobid

# execute code
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1
time ./cov_omp > output/$PBS_JOBID.out
