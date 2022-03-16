/**
 * Modern Computing Technologies, WUT, 2022
 * 
 * Author: Jakub Zielinski
 * Date: 16.03.2022
 * 
 * Compilation command:
 *      gcc covariance.c -o cov_mpi -lm 
 * */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "mct_utils.h"

#define NELEMENTS 134217728

int main( int argc , char ** argv ) 
{
    int i, j;
    
    // allocate memory
    double *var[5];
    for(i=0; i<5; i++) cppmallocl(var[i],NELEMENTS,double);     

    
    // Read binary data
    b_t(); // start timing
    FILE * fp;
    char file_name[128];
    for(i=0; i<5; i++)
    {
        sprintf(file_name, "/home2/archive/MST-2022/lab1/var%d.dat", i+1);
        printf("# READING DATA: `%s`...\n", file_name);
        
        fp = fopen(file_name, "rb");
        if(fp==NULL) { printf("# ERROR: Cannot open `%s`!\n", file_name); return EXIT_FAILURE; }
        
        size_t ftest = fread(var[i], sizeof(double), NELEMENTS, fp);
        if(ftest!=NELEMENTS) { printf("# ERROR: Cannot read `%s`!\n", file_name); return EXIT_FAILURE; }
        
        fclose(fp);
    }
    double tio = e_t(); // stop timing
    printf("# READ TIME: %f sec\n", tio);


    long int idx[10] = {10,13421673,25501328,41606496,53677091,73818750,83214911,93952210,106032910,132875451};
    for(i=0; i<5; i++) 
      for(j=0; j<10; j++)
	printf("var[%d][%10ld]=%16.8g\n", i,idx[j],var[i][idx[j]]);
    
    /* 
    // generate additional random variables v_6, ..., v_10 and compute covariance matrix
    // make computation only for elements cov(i,j) where i<=j (see printing statment)
    double cov[10][10]; // storage for covaraince matrix
    b_t(); // start timing
    
    // TODO: add your code here!
    
    double tcmp = e_t(); // stop timing
    printf("# COMPUTATION TIME: %f sec\n", tcmp);
    
    // print results
    for(i=0; i<10; i++) for(j=0; j<=i; j++)
    {
        printf("cov(%2d,%2d)=%16.8g\n", i+1, j+1, cov[i][j]);
    }
    */
    
    return( EXIT_SUCCESS ) ;
}
