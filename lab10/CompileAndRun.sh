#!/bin/bash
SOURCE="laplace2d-cuda.cu"
FLAGS="-lcufft -lcudart -lm" #"-I/usr/local/cuda/inc -L/usr/local/cuda/lib"

# module load openmpi-gcc447-Cuda90/3.1.1

while test $# -gt 0; do
    case "$1" in
	-h|--help)
	    echo "Script for compiling and runnig code for Task 5 (MCT)"
	    echo "Available options:"
	    echo "-h, --help        prints brief description of arguments"
	    break
	    ;;
	*)
	    break
	    ;;
    esac
done

rm -f laplace2d-cuda
# gcc -std=c99 $SOURCE -o laplace2d $FLAGS
nvcc -arch sm_35 -O3 $SOURCE -o laplace2d-cuda $FLAGS

if [ -e laplace2d-cuda ]
then
    ./laplace2d-cuda > output.txt
else
    echo "# Compilation failed!"
fi
