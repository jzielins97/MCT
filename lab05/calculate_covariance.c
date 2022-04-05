/**
 * Modern Computing Technologies, WUT, 2022
 * 
 * Author: Jakub Zieliński
 * Date: 30.03.2022
 * 
 * Compilation command:
 *      gcc -shared -o libcov.so -fPIC calculate_covariance.c
 *
 *
 * Covariance matrix calculation in functions
 * to embed in Python
 * 
 * 
 * */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h> // include openMP library

#define CHUNKSIZE 700000

void read_data_test(long int N, int id,  double* var){
  printf("\t%f\n",var[id*N+10]); // printing 10th element
  /* for(int in=0; in<N; in++){ */
  /*   if(in%10000 == 0){ */
  /*     printf("\t%f\n",var[id*N+in]); */
  /*   } */
  /* } */
}

// function for taking 10 * 1D array because 2D does not work, yay
void calculate_other_variables(long int N,
			       double* var0,
			       double* var1,
			       double* var2,
			       double* var3,
			       double* var4,
			       double* var5,
			       double* var6,
			       double* var7,
			       double* var8,
			       double* var9){
  printf("Inside calculate_covariance function in c (double** var)\n");
  long int in=0;
  int chunk = CHUNKSIZE;
  int nthreads = 1;
  // get number of threads
#pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  chunk = ceil(1.0 * N / nthreads ); //CHUNKSIZE;
  printf("chunks are of size %d (num. of threads = %d)\n",chunk, nthreads);
  printf("10th elements are:\n");
  printf("\t%f",var0[10]);
  printf("\t%f",var1[0]);
  printf("\t%f",var2[10]);
  printf("\t%f",var3[10]);
  printf("\t%f\n",var4[10]);

  printf("Starting for loop\n");
  register double var0i, var1i, var2i, var3i, var4i, var5i, var6i, var7i, var8i, var9i;  
  #pragma omp for schedule(dynamic,chunk)
  for(in = 0; in < N; in++){
    // read - 5
    var0i = var0[in];
    var1i = var1[in];
    var2i = var2[in];
    var3i = var3[in];
    var4i = var4[in];

    // compute
    var5i = sin(var1i) + sin(var0i);
    var6i = exp(var2i) - exp(-1.*var4i);
    var7i = sin(var3i)*cos(var0i) + cos(var3i)*sin(var2i);
    var8i = hypot(var2i, var1i);
    var9i = cbrt(var3i);

    // write - 5
    var5[in] = var5i;
    var6[in] = var6i;
    var7[in] = var7i;
    var8[in] = var8i;
    var9[in] = var9i;
    
  }// end of for loop
  printf("Finished for loop\n");
}

void calculate_covariance(long int N,  // number of elements
			  double* var, // pointer to 1D array (10*N) with variables
			  double* avg,
			  double* cov) // pointer to 1D array (10) with averages
{
  int i, j, in;
  int chunk = CHUNKSIZE;
  int nthreads = 1;
  double priv_cov;
  #pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  for(i=0; i<10; i++){
    for(j=0; j<10; j++){
      priv_cov = 0;
      #pragma omp for schedule(dynamic,chunk)
      for(in=0; in<N; in++){
	priv_cov += (var[i*N+in] - avg[i]) * (var[j*N+in] - avg[j]);
      }
      #pragma omp critical
      cov[i*10 + j] += priv_cov / (N-1);
    }
  }
}

// function for taking 1D array and iterating over it as 2D
void calculate_covariance_1Darray(long int N, double* _var, double* _cov){
  printf("Inside calculate_covariance function in c (double** var)\n");
  long int i, j, in=0;
  // create pointers to specific memory in the array
  double* var[10];
  double* cov[10];
  for(int i=0; i<10; i++){
    var[i] = &_var[i*N];
    cov[i] = &_cov[i*10];
  }
  double avg_var[10] = { 0 };
  int chunk = CHUNKSIZE;
  int nthreads = 1;
  // get number of threads
  #pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  chunk = ceil(1.0 * N / nthreads ); //CHUNKSIZE;
  printf("chunks are of size %d (num. of threads = %d)\n",chunk, nthreads);
  printf("10th elements are:\n");
  printf("\t%f",var[0][10]);
  printf("\t%f",var[1][10]);
  printf("\t%f",var[2][10]);
  printf("\t%f",var[3][10]);
  printf("\t%f\n",var[4][10]);

  printf("Starting for loop\n");
  double avg_priv[10] = { 0 };
  register double var0i, var1i, var2i, var3i, var4i, var5i, var6i, var7i, var8i, var9i;  
  #pragma omp for schedule(dynamic,chunk)
  for(in = 0; in < N; in++){
    // read - 5
    var0i = var[0][in]; // var[0*N+in];
    var1i = var[1][in]; // var[1*N+in];
    var2i = var[2][in]; // var[2*N+in];
    var3i = var[3][in]; // var[3*N+in];
    var4i = var[4][in]; // var[4*N+in];
    
    // write - 5
    avg_priv[0] += var0i;
    avg_priv[1] += var1i;
    avg_priv[2] += var2i;
    avg_priv[3] += var3i;
    avg_priv[4] += var4i;

    // compute
    var5i = sin(var1i) + sin(var0i);
    var6i = exp(var2i) - exp(-1.*var4i);
    var7i = sin(var3i)*cos(var0i) + cos(var3i)*sin(var2i);
    var8i = hypot(var2i, var1i);
    var9i = cbrt(var3i);

    // write - 5
    var[5][in] = var5i;
    var[6][in] = var6i;
    var[7][in] = var7i;
    var[8][in] = var8i;
    var[9][in] = var9i;
    
    // write - 5
    avg_priv[5] += var5i;
    avg_priv[6] += var6i;
    avg_priv[7] += var7i;
    avg_priv[8] += var8i;
    avg_priv[9] += var9i;
    if(in%10000 == 0) printf("%d\n",i);
  }// end of for loop
  printf("Finished for loop\n");
  #pragma omp critical
  for(i = 0; i<10; i++) avg_var[i] += avg_priv[i] / N;

  printf("Calculating covariance\n");
  double priv_cov = 0;
  for(i=0; i < 10; i++){
    for(j = 0; j<=i; j++){
      priv_cov=0;
      #pragma omp for schedule(dynamic,chunk)
      for(in=0; in<N; in++)
	priv_cov += (var[i][in] - avg_var[i]) * (var[j][in] - avg_var[j]);
      #pragma omp critical
      cov[i][j] += priv_cov / (N - 1);
    }
  }
  printf("finished calculations\n");
}

  
    
    
