#!/bin/bash
SOURCE="laplace2d.c"
FLAGS="-lfftw3 -lm"

# module load openmpi-gcc447-Cuda90/3.1.1

while test $# -gt 0; do
    case "$1" in
	-h|--help)
	    echo "Script for compiling and runnig code for Task 5 (MCT)"
	    echo "Available options:"
	    echo "-h, --help        prints brief description of arguments"
	    echo "-omp              flag for compiling code with OpenMP multithreading. Standard option is using a serial code"
	    echo "-nt               specify value for OMP_NUM_THREADS"
	    exit 0
	    ;;
	-omp)
	    shift # removes first argument from the $#
	    echo "Compiling with OpenMP"
	    SOURCE="laplace2d-omp.c"
	    FLAGS="-lfftw3_threads -lfftw3 -lm -lpthread -fopenmp -lgomp"
	    ;;
	-nt)
	    shift
	    if [ $# -gt 0 ]; then 
		export OMP_NUM_THREADS=$1
		shift
	    else
		export OMP_NUM_THREADS=1
		echo "No number of threads specified, using 1"
	    fi
	    ;;
	*)
	    break
	    ;;
    esac
done

rm -f laplace2d
gcc -std=c99 $SOURCE -o laplace2d $FLAGS

if [ -e laplace2d ]
then
    time ./laplace2d
else
    echo "# Compilation failed!"
fi
