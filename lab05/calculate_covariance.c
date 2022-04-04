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

void read_data_test(long int N, double* var){
  printf("Inside read_data_test function in c (double* var)\n");
  for(int in=0; in<N; in++){
    if(in%N == 10){
      printf("\t%f\n",var[in]);
      break;
    }
  }
}

void calculate_covariance(long int N, double** var, double** cov){
  printf("Inside calculate_covariance function in c (double** var)\n");
  double avg_var[10] = { 0 };
  long int i, j, in=0;
  int chunk = CHUNKSIZE;
  int nthreads = 1;
  printf("var: ");
  do{
    if(in%N == 0){
      printf("\n***************\n");
      printf("(%d, %f), ", in, var[in]);
    }    
    if(in%N == 10) printf("(%d, %f), ", in, var[in]);
    in++;
  }while(var[in] != NULL);
  printf("\nTotal number of elements: %d\n", in);
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
  printf("\t%f",var[4][10]);

#pragma omp parallel shared(var) private(in,i,j)
  {
    printf("Starting for loop");
    double avg_priv[10] = { 0 };
    register double var0i, var1i, var2i, var3i, var4i, var5i, var6i, var7i, var8i, var9i;  
    #pragma omp for schedule(dynamic,chunk)
    for(in = 0; in < N; in++){
      // read - 5
      var0i = var[0][in];
      var1i = var[1][in];
      var2i = var[2][in];
      var3i = var[3][in];
      var4i = var[4][in];
    
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
    }// end of for loop
    printf("Finished for loop");
    #pragma omp critical
    for(i = 0; i<10; i++) avg_var[i] += avg_priv[i] / N;
  }//end of openMP

  double priv_cov = 0;
#pragma omp parallel private(in) reduction(+:priv_cov)
  {
    for(i=0; i < 10; i++){
      for(j = 0; j<=i; j++){
	priv_cov=0;
        #pragma omp for schedule(dynamic,chunk)
	for(in=0; in<N; in++)
	  priv_cov += (var[i][in] - avg_var[i]) * (var[j][in] - avg_var[j]);
	cov[i][j] = priv_cov / (N - 1);
      }
    }
  }//end of openMP
  printf("finished calculations\n");
}


  
    
    
