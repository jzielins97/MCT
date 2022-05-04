#!/bin/bash

if [ "$LD_LIBRARY_PATH" != "/usr/local/lapack-3.9.0-gcc721/lib64/" ]; then
    echo "Exporting library path for lapacke"
    export LD_LIBRARY_PATH=/usr/local/lapack-3.9.0-gcc721/lib64/
fi

gcc eigenvalue.c -o eigenvalue -O3 -lcblas -llapacke -I/usr/local/lapack-3.9.0-gcc721/include/ -I/usr/include/ -L/usr/lib64/atlas/ -L/usr/local/lapack-3.9.0-gcc721/lib64/ -lm

if test -e eigenvalue
then
    ./eigenvalue
fi