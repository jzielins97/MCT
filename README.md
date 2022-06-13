# MCT
Modern Computing Technologies (2022)


## Lab 1 (Task 1):

Modern Computing Technologies, WUT, 2022

Author: [YOUR_NAME]
Date: [DATE_OF_SENDIG_THE_REPORT]

Compilation command:
gcc covariance.c -o covariance -lm

TASK DESCRIPTION

Write a code that computes covariance matrix of random variables v_1, ..., v_10.
Compute covariance matrix using formula:
https://en.wikipedia.org/wiki/Covariance#Calculating_the_sample_covariance

1. Random variables v_1, ..., v_5 are provided in binary files var1.dat, ..., var5.dat respectively.
Files are located in path:
/home2/archive/MCT-2022/lab1
Each series (file) contain 134217728 elements (doubles).

2. Generate random variables v_6, ..., v_10 according formulas:\
v<sub>6</sub> = sin(v<sub>2</sub>) + sin(v<sub>1</sub>);\
v<sub>7</sub> = exp(v<sub>3</sub>) - exp(-v<sub>5</sub>);\
v<sub>8</sub> = sin(v<sub>4</sub>)*cos(v<sub>1</sub>) + cos(v<sub>4</sub>)*sin(v<sub>3</sub>);\
v<sub>9</sub> = hypot(v<sub>3</sub>, v<sub>2</sub>);\
v<sub>10</sub>= cbrt(v<sub>4</sub>);\

3. Measure time needed for loading data and for computation of the covariance matrix,
for example:
 READ TIME: 34.508760 sec
...
 COMPUTATION TIME: 105.915668 sec

4. As output provide values of covariance matrix in form (put them to the report):
cov( 1, 1)= 0.33330407
cov( 2, 1)= 0.33330627
cov( 2, 2)= 0.354141
.... ....
cov(10, 9)= 0.30952735
cov(10,10)= 0.14192182

5. To speed up you work start with provided template below.

6. Deadline for submission of report is: 09-03-2022


## Lab 2 (Task 2):


Speed up your Task1 code for the covariance matrix calculation using openMP.

Execute calculations for following number of threads:
        1, 2, 4, 8, 16, 20, 30, 40.

In the report provide simple plot showing the scaling.

EDIT: USE job.sh SCRIPT TO CALL A JOB. DO NOT USE NODE 72 IN INTERACTIVE MODE.

## Lab 3-4 (Task 3):
### PART I
Start from the code: /home2/archive/MCT-2022/lab3

Convert into MPI form part that is responsible for reading data (arrays: var1-var5). Code should read data from files and place it into operating memory in block distribution (use: getNip(), and get_i0() from attached mct_utils.h in order to determine chunk sizes).

Use MPI_Send() and MPI_Recv() routines. Check correctness of data (by checking values of selected elements of array, line: 52 of the template) when code is executed on full dwarf system. For this use template of submission script. 

Do not write report at this stage – in next lab you will add remaining part of the code that computes the covariance matrix.

NOTE: The proper way of reading data in this case is to use MPI I/O, which we will discuss later. Aim if this lab is to learn how to operate with MPI_Send( … ) and MPI_Recv( … ) routines

### PART II
​Start from the code you have prepared for Part I

Using MPI collective routines speed-up the code for covariance computation  MPI_Allreduce().

Execute code with: 1, 2, 4, 8, 16, 20, 30, 40, ... processes and perform time measurements (mimum on 2 nodes).

On a single plot compare timing obtained for openMP code and MPI code (data for openMP code are limited only only to interval 1-40).

Add to the plot line corresponding to the ideal scaling:
t_n = t1 / n,
 where: t_n – time measured using n computing units, 
t1 – time measured for serial code (single computing unit) 
n – number of computing units
Check correctness of computation – compare results with reference code (serial) .

## Lab 5 (Task 4):
Rewrite your code for covariance matrix (OpenMP version) to Python + C  framework.

Identify computationally intensive parts and implement them in C  language and use Python to create 'user friendly' code/API.

---------------------------------------------
REPORT:
    C and python code + terminal output

## Lab 6 (Task 5)
Starting from provided template write a code that computes laplace (in 2D) of function f(x,y):\
    <img src="https://latex.codecogs.com/svg.image?\bg{white}F(x,y)&space;=&space;\Delta_{2D}&space;f(x,y)&space;=&space;\frac{\partial^2&space;f}{\partial&space;x^2}&space;&plus;&space;\frac{\partial^2&space;f}{\partial&space;y^2}" title="https://latex.codecogs.com/svg.image?\large \bg{white} F(x,y) = \Delta_{2D} f(x,y) = \frac{\partial^2 f}{\partial x^2} + \frac{\partial^2 f}{\partial y^2}" />

Test you code using function (already in the template: /home2/archive/MCT-2022/lab6 )\
    <img src="https://latex.codecogs.com/svg.image?\large&space;\bg{white}f(x,y)&space;=&space;exp(Ax^2&space;&plus;&space;By^2&space;&plus;&space;Cxy)" title="https://latex.codecogs.com/svg.image?\large \bg {white}f(x,y) = exp(Ax^2 + By^2 + Cxy)" />

Use openMP technology (see section http://www.fftw.org/fftw3_doc/Multi_002dthreaded-FFTW.html) 

Execute runs with openMP for nx=ny=4112

---------------------------------------------
REPORT: Code + Scaling plot (processses: 1, 2, 4, 8, 16, 20) + output of the test runs
Weight for this assignment is 1.5

---------------------------------------------------------
TIPS:
Compute forward FT of f(x,y) to obtain f(k_x, k_y).
From that Construct F(k_x, k_y) = (-k_x^2 - k_y^2) f(k_x, k_y)/(n_x n_y).
Then compute backward FT of F(k_x, k_y) to obtain F(x,y) = \lap_{2D} f(x,y)

##Lab 7 (Task 6):
Using provided template in /home2/archive/MCT-2022/lab7 write a code that computes all eigenvalues and eigenvectors of a given matrix_H() function.

Check correctness of computation by comparing with prediction for harmonic oscilator: 
https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator
E=\hbar \omega (n+ 1/2)

Plot first 5 eigenvectors, which qualitatively (to the sign) should agree with eigenfunctions with:
https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator#/media/File:HarmOsziFunktionen.png

-----------------------
REPORT:
code
comparison of computed eigenvalues with analytic solution
plots of 5 eigenvecotrs 
estimate of what is the maximum matrix size that you can diagonalize in tme shorter than 10 min

-----------------------
TIPS:
Note that in line 25 of the given template one can find a function that should be used in order to diagonalize the Hamiltoniam Matrix.

## Lab 8-9 (Task 7):
Start with the code that you constructed during lab 3 and rewrite the part of the code that is responsibe for reading the files (var1-var5) using MPI I/O functions, that were presented during the lecture. 

Compare timings with the previous solution (MPI_Send/Recv and MPI I/O).
Additionally add block of code that will write into files var6-var10, again using MPI I/O functions.

-----------------------
REPORT:
+code
+output of the code (with covariance matrix)
+comparison of reading time of previous solution with the MPI I/O.
+write time of var6-var10
+plots showing measured bandwidth for read and write (MB/sec) as a function of number of MPI processes

## Lab 10-11 (Task 8)
Instructions
Write a code that computes laplace (in 2D) of function f(x,y) using cuFFT and CUDA
Start with the code constructed for lab 5. 

Use following links as a hints:
https://docs.nvidia.com/cuda/cufft/index.html
https://docs.nvidia.com/cuda/cufft/index.html#fftw-conversion-guide
https://docs.nvidia.com/cuda/cufft/index.html#cufft-code-examples
Execute test runs with cuda for nx=ny=4112

-------------------------
REPORT: 
+code
+output
+timngs for computation and data transfers

## Lab 12-13 (Task 9)
Instructions
Preparation for lab 13:
Register to PL-Grid: https://portal.plgrid.pl/registration/form
For more info see: https://docs.cyfronet.pl/display/PLGDoc/User+Manual#UserManual-Basicinformation
To confirm affilition, you can use my data as “Supervisor’s affiliation”:
Name: Gabriel
Surname: Wlazlowski
E-mail: gabriel.wlazlowski@pw.edu.pl
OPI: 225935

    2. Create trial grant
    3. Activate service: access to Zeus supercomputer, or access to Eagle supercomputr, or other
    4. Try to login to the granted supercomputer


Task for lab 13:
    1. Execute strong scaling of MPI code (you can use inversion-scalapack.c)
    2. Execute strong scaling using one of computers provided by PL-Grid (via trial grant)
    3. Suppose that your problem requires 100 inversions of matrix of size 25,000 x 25,000. For such problem propose parameters of runs (how many nodes/cores), estimate time-to-solution, and estimate needed budget (CPU-hours),


-----------------------
In the report provide:
→ compilation script / command and sample of your submission script,
→ strong scaling plot,
→ Proposition of execution parameters with time-to-solution and cost estimation with short justification, why you selected these parameters.


Time-to-solution: time needed by a computer to solve you problem (hr, min, sec, …) strongly depends on running parameters like number of processes. Computing cost = Time-to-solution x number of processes (CPU-hr, Node-hr, GPU-hr…)

Example:
if your code was running for 1 hour using 100 CPUs then computing cost was 100 CPU-hr.
