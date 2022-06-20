#!/bin/bash
##dodanie informacji o grancie obliczeniowym:
## Name of the grant
#SBATCH -A plgjzielins972022a
## Name of the job
#SBATCH -J jzielins97-scal
## Resources: number of nodes, number of processes per node
#SBATCH -N 10
#SBATCH --ntasks-per-node=24
## Maximum memory for one process (default 5GB)
#SBATCH --mem-per-cpu=3GB
## Maximum time for execution
#SBATCH --time=1:00:00
## Queue
#SBATCH -p plgrid
## Standard output/error
#SBATCH --output="./slurm/%j.out"
##SBATCH --error="error.err"

srun /bin/hostname

#script arguments
nc=1 # no compilation flag
run=1 # flag to run program
#-------- Ladowanie modules -------------------->
## Load module OpenMPI
#module load plgrid/tools/openmpi/2.1.1-gcc-6.4.0
#module add plgrid/tools/openmpi/4.1.2-nvhpc-22.3
## Load modulu scalapack
module add plgrid/libs/scalapack/2.2.0-nvhpc-22.3
## Load module atlas
module add plgrid/libs/atlas/3.10.3
## Load module OpenBlas
#module add plgrid/libs/openblas/0.3.15
#-----------------------------------------------<


#--------- Set compilation variables ----------->
OPENMPI_LIB="/net/software/local/openmpi/2.1.1-gcc-6.4.0/lib"
#SCALAPACK_LIB="/net/software/local/scalapack/2.2.0-nvhpc-22.3/lib"
ATLAS_LIB="/net/software/local/atlas/3.10.3"
OPENBLAS_LIB="/net/software/local/openblas/0.3.15/lib"
FLAGS="-L$OPENMPI_LIB -lscalapack"
SOURCE="inversion-scalapack.c"
TARGET="inversion-scalapack"
NP=240
#-----------------------------------------------<


# go to right directory
cd $SLURM_SUBMIT_DIR
pwd

#------- COMPILATION --------------------------->
if [ $nc -eq 0 ]
then
	# compile code
	if [ -e $SOURCE ]
	then
		rm -f $TARGET
		mpicc $SOURCE -o $TARGET $FLAGS
		if [ ! -e $TARGET ]
		then
			echo "#ERROR: compilation failed!"
			exit -1
		else
			echo "Compilation finished"
		fi
	else
		echo "#ERROR: there is no $SOURCE in $(pwd)"
		exit -1
	fi
fi
#-----------------------------------------------<

#-------- run program with MPI ----------------->
if [ $run -eq 1 ]
then
	if [ -e $TARGET ]
	then
		mpirun ./$TARGET > OUTPUTS/$NP.txt
	else
		echo "#ERROR: there is no target $TARGET in $(pwd)"
	fi
fi
#----------------------------------------------<
############## END SCRIPT ######################

# In order the s script use:
# sbatch job.sh

# In order to check status of job use:
# squeue

# In order to delete job use:
# scancel jobid
