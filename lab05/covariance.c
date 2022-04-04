/**
 * Modern Computing Technologies, WUT, 2022
 * 
 * Author: Jakub Zieliński
 * Date: 09.03.2022
 * 
 * Compilation command:
 *      export OMP_NUM_THREADS=20
 *      gcc cov_omp.c -O3 -o cov_omp -lm -fopenmp
 *
 *
 * TASK DESCRIPTION
 *   Speed up the covariance code using openMP.
 *   In the report provide a graph of time vs number of threads
 *
 * Deadline for submission of report is: 16-03-2022
 * */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h> // include openMP library

#include "mct_utils.h"


#define NELEMENTS 134217728

int main(int argc, char **argv)
{
  FILE* fout = fopen("time.txt","a"); // file for storring timew
  int nthreads; // number of threads in OpenMP
  if(fout == NULL) {
    printf("# Error: cannot open the file");
    return EXIT_FAILURE;
  }
  int i, j, in;

    // allocate memory
    double *var[10];
    for (i = 0; i < 10; i++)
        cppmallocl(var[i], NELEMENTS, double);

    // Read binary data
    b_t(); // start timing
    FILE *fp;
    char file_name[128];
    for (i = 0; i < 5; i++)
    {
        sprintf(file_name, "/home2/archive/MCT-2021/lab1/var%d.dat", i + 1);
        printf("# READING DATA: `%s`...\n", file_name);

        fp = fopen(file_name, "rb");
        if (fp == NULL)
        {
            printf("# ERROR: Cannot open `%s`!\n", file_name);
            return EXIT_FAILURE;
        }

        size_t ftest = fread(var[i], sizeof(double), NELEMENTS, fp);
        if (ftest != NELEMENTS)
        {
            printf("# ERROR: Cannot read `%s`!\n", file_name);
            return EXIT_FAILURE;
        }

        fclose(fp);
    }
    double tio = e_t(); // stop timing
    printf("# READ TIME: %f sec\n", tio);

    // get number of threads
#pragma omp parallel
    {
      nthreads = omp_get_num_threads();
    }

    // generate additional random variables v_6, ..., v_10 and compute covariance matrix
    // make computation only for elements cov(i,j) where i<=j (see printing statment)
    double cov[10][10] = { 0 }; // storage for covaraince matrix
    double sum_var[10] = { 0 };
    double avg_var[10] = { 0 };
    double var0i, var1i, var2i, var3i, var4i, var5i, var6i, var7i, var8i, var9i; //memory optimalization
    b_t();              // start timing

    //starting openMP
    int chunk = ceil(1.0 * NELEMENTS / nthreads); // how many elements to give to each thread
#pragma omp parallel shared(var) private(in,i,j)
    {
    //generate additional random variables
      double priv_avg[10] = { 0 };
      #pragma omp for schedule(dynamic,chunk)
      for(in = 0; in < NELEMENTS; in++){
	/* priv_sum[0] += var[0][in]; */
	/* priv_sum[1] += var[1][in]; */
	/* priv_sum[2] += var[2][in]; */
	/* priv_sum[3] += var[3][in]; */
	/* priv_sum[4] += var[4][in]; */

	/* var[5][in] = sin(var[1][in]) + sin(var[0][in]); */
	/* var[6][in] = exp(var[2][in]) - exp(-1.*var[4][in]); */
	/* var[7][in] = sin(var[3][in])*cos(var[0][in]) + cos(var[3][in])*sin(var[2][in]); */
	/* var[8][in] = hypot(var[2][in], var[1][in]); */
	/* var[9][in] = cbrt(var[3][in]); */

	/* priv_sum[5] += var[5][in]; */
	/* priv_sum[6] += var[6][in]; */
	/* priv_sum[7] += var[7][in]; */
	/* priv_sum[8] += var[8][in]; */
	/* priv_sum[9] += var[9][in]; */

	// read - 5
	var0i = var[0][in];
	var1i = var[1][in];
	var2i = var[2][in];
	var3i = var[3][in];
	var4i = var[4][in];

	// write - 5
	priv_avg[0] += var0i;
	priv_avg[1] += var1i;
	priv_avg[2] += var2i;
	priv_avg[3] += var3i;
	priv_avg[4] += var4i;

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
	priv_avg[5] += var5i;
	priv_avg[6] += var6i;
	priv_avg[7] += var7i;
	priv_avg[8] += var8i;
	priv_avg[9] += var9i;
  
      }
      #pragma omp critical
      for(i = 0; i<10; i++) avg_var[i] += priv_avg[i];
    }// end of openMP

    for(i = 0; i<10; i++) avg_var[i] = avg_var[i] / NELEMENTS;

    // calculate covariance
    double priv_cov = 0;
#pragma omp parallel private(in) reduction(+:priv_cov)
    {

      for (i = 0; i < 10; i++){
	for(j = 0; j <= i; j++){
	  priv_cov = 0;
          #pragma omp for schedule(dynamic,chunk)
	    for( in=0; in<NELEMENTS; in++)
	      priv_cov += (var[i][in] - avg_var[i]) * (var[j][in] - avg_var[j]);
	  cov[i][j] = priv_cov / (NELEMENTS - 1);
	}
      }
    }//end of openMP
	

    double tcmp = e_t(); // stop timing
    printf("# COMPUTATION TIME: %f sec\n", tcmp);
    printf("# number of threads = %d\n", nthreads);
    fprintf(fout,"%d\t%lf\n",nthreads,tcmp);
    fclose(fout);

    // print results
    for (i = 0; i < 10; i++)
        for (j = 0; j <= i; j++)
            printf("cov(%2d,%2d)=%16.8g\n", i + 1, j + 1, cov[i][j]);

    return (EXIT_SUCCESS);
}
