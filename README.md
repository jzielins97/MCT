# MCT
Modern Computing Technologies (2022)


## Lab 1:

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

2. Generate random variables v_6, ..., v_10 according formulas:
v_6 = sin(v_2) + sin(v_1);
v_7 = exp(v_3) - exp(-1.*v_5);
v_8 = sin(v_4)*cos(v_1) + cos(v_4)*sin(v_3);
v_9 = hypot(v_3, v_2);
v_10= cbrt(v_4);

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


## Lab 2:


Speed up your Task1 code for the covariance matrix calculation using openMP.

Execute calculations for following number of threads:
        1, 2, 4, 8, 16, 20, 30, 40.

In the report provide simple plot showing the scaling.

EDIT: USE job.sh SCRIPT TO CALL A JOB. DO NOT USE NODE 72 IN INTERACTIVE MODE.

## Lab 3 & 4:
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



