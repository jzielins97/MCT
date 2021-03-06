#!/bin/bash
SOURCE="covariance.c"
TARGET="covariance_mpiio"
FLAGS="-O3 -lm"
NP=1 # number of processes with MPI
fNC=0  # flag for compiling or not
fLAB3=0 #flag to identify which program is running

while test $# -gt 0; do
    case "$1" in
	-h|--help)
	    echo "Script for compiling and runnig code for Task 5 (MCT)"
	    echo "Available options:"
	    echo "-h, --help        prints brief description of arguments"
	    echo "-np               specify number of processes"
	    echo "-nc               run without recompiling"
	    echo "-lab3             run program from Lab3 (send/recv routine with MPI)"
	    exit 0
	    ;;
	-nc)
	    shift
	    fNC=1
	    ;;
	-np)
	    shift
	    if [ $# -gt 0 ]; then 
		NP=$1
		shift
	    else
		echo "No number of threads specified, using 1"
	    fi
	    ;;
	-lab3)
	    shift
	    echo "Running program from the Lab3 (send/recv routine)"
	    fLAB3=1
	    SOURCE="covariance_lab3.c"
	    TARGET="covariance_mpi"
	    ;;
	*)
	    break
	    ;;
    esac
done

if [ $fNC -eq 0 ]
then
    echo "Compiling code..."
    rm -f $TARGET
    mpicc $SOURCE -o $TARGET $FLAGS
    if [ ! -e $TARGET ]
    then
	echo "# Compilation failed!"
    else
	echo "finished"
    fi
fi


if [ -e $TARGET ]
then
    if [ $fLAB3 -eq 1 ]
    then
	mpirun -np $NP ./$TARGET > OUTPUTS/lab3/output_$NP.txt
    else
	mpirun -np $NP ./$TARGET > OUTPUTS/output_$NP.txt
    fi
else
    echo "# No file called $TARGET"
fi
