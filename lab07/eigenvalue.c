/**
 * Task 6 template
 * 
 * gcc eigenvalue.c -o eigenvalue -O3 -I/usr/include/ -L/usr/lib64/atlas/ -l cblas -llapack -lm
 * */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "mct_utils.h"

// use C interface
// convenient documentation is provided here:
// https://software.intel.com/en-us/mkl-developer-reference-c
#include <cblas.h>

// Invoke directly Fortran functions
// see: https://software.intel.com/en-us/mkl-developer-reference-fortran

// call dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)

#define N 400
/**
 * Function returns value of matrix element A_{ij}
 * i,j iterates from 1 to M as in standard mathematical notation
 * */
double matrix_H(int k, int l)
{
#define omega 0.001
#define a 1.0
#define hbar 1.0
#define m 1.0

    double K, V;
    if (k == l)
    {
        K = pow(hbar * M_PI / a, 2) / (6. * m) * (1.0 + 2.0 / (N * N));
        V = 0.5 * m * pow(omega * a * (k - N / 2), 2);
    }
    else
    {
        K = pow(hbar * M_PI / a, 2) / (1. * m * N * N) * pow(-1, k - l) / pow(sin(M_PI * (k - l) / N), 2);
        V = 0.0;
    }

    return K + V;

#undef a
#undef m
}

#define IDX1(i, j) (j - 1) * N + (i - 1)

int main()
{

    double *H; // matrix
    cppmallocl(H, N * N, double);

    int i, j;
    double rt;

    // set matrix
    for (i = 1; i <= N; i++)
        for (j = 1; j <= N; j++)
            H[IDX1(i, j)] = matrix_H(i, j);

    // Compute eigen values and eigen vectors
    b_t();

    // TODO

    double rt = e_t();
    printf("Computation time: %fsec\n", rt);

    return 1;
}
