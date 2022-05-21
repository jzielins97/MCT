#!/bin/bash

for np in 40 40 20 20 16 16 8 8 4 4 2 2 1 1
do
    echo "NP=$np"
    echo "    -Running code with mpiio"
    ./CompileAndRun.sh -nc -np $np
    echo "    -Running code with send/recv"
    ./CompileAndRun.sh -nc -np $np -lab3
done

