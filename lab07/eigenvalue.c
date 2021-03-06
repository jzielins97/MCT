/**
 * Task 6 template
 * 
 * export LD_LIBRARY_PATH=/usr/local/lapack-3.9.0-gcc721/lib64/
 *
 * gcc eigenvalue.c -o eigenvalue -O3 -lcblas -llapacke -I/usr/local/lapack-3.9.0-gcc721/include/ -I/usr/include/ -L/usr/lib64/atlas/ -L/usr/local/lapack-3.9.0-gcc721/lib64/ -lm
 *
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
#include <lapacke.h>

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
  /* int N = 400; */
  double *H; // matrix
  cppmallocl(H, N * N, double);
  
  int i, j;
  double rt;
  
  // set matrix
  for (i = 1; i <= N; i++)
    for (j = 1; j <= N; j++)
      H[IDX1(i, j)] = matrix_H(i, j);
  
  double* W; // vector of eigen values
  cppmallocl(W, N, double);
  
  // using subroutine for calculating eigenvalues of a symetric matrix:
  // SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO )
  // http://www.netlib.org/lapack/double/dsyevd.f
  char jobz = 'V'; // V for storing eigenvectors, N to ignore
  char uplo = 'L'; // storing upper (U) or lower (L) triangle of A
  // leading dimensions:
  // it is a number of elements in a row(if using row-major notaion)
  int ldh = N; // using NxN matrix
  int info; // output of the subroutine (0 - finished succesfuly)
  
  // Compute eigen values and eigen vectors
  b_t();
  
  // TODO
  info = LAPACKE_dsyevd(CblasRowMajor, jobz, uplo, N, H, ldh, W);
  // LAPACKE_dsyevd (CblasRowMajor, jobz, uplo, N, H, ldh, W, WORK,  ldw, IWORK, ldi);
  if(info!=0) { printf("Error: LAPACKE_dsyevd=%d\n", info); return 1;}
  
  rt = e_t();
  printf("Computation time: %fsec\n", rt);

  FILE* fValues = fopen("eigenvalues.txt","w"); // file for storring eigenvalues
  FILE* fVectors = fopen("eigenvectors.txt","w"); // file for storring vectors
  fprintf(fVectors,"%d\n",N); // first value in the file is dimention of the vectors
  for(int i=0; i<N; i++){
    fprintf(fValues,"%lf\n",W[i]);
    for(int j=0; j<ldh;j++)
      fprintf(fVectors,"%g ", H[j*N+i]);
    fprintf(fVectors,"\n");
  }
  fclose(fValues);
  fclose(fVectors);
  return 1;
}

