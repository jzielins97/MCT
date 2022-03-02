/**
 * Modern Computing Technologies, WUT, 2022
 * 
 * Author: [YOUR_NAME]
 * Date: [DATE_OF_SENDIG_THE_REPORT]
 * 
 * Compilation command:
 *      gcc covariance.c -o covariance -lm 
 *
 *
 * TASK DESCRIPTION
 * 
 * Write a code that computes covariance matrix of random variables v_1, ..., v_10.
 * Compute covariance matrix using formula:
 *      https://en.wikipedia.org/wiki/Covariance#Calculating_the_sample_covariance
 * 
 * 1. Random variables v_1, ..., v_5 are provided in binary files var1.dat, ..., var5.dat respectively.
 *    Files are located in path:
 *        /home2/archive/MCT-2022/lab1
 *    Each series (file) contain 134217728 elements (doubles).
 * 
 * 2. Generate random variables v_6, ..., v_10 according formulas:
 *       v_6 = sin(v_2) + sin(v_1);
 *       v_7 = exp(v_3) - exp(-1.*v_5);
 *       v_8 = sin(v_4)*cos(v_1) + cos(v_4)*sin(v_3);
 *       v_9 = hypot(v_3, v_2);
 *       v_10= cbrt(v_4);
 * 
 * 3. Measure time needed for loading data and for computation of the covariance matrix,
 *    for example:
 *    # READ TIME: 34.508760 sec
 *    ...
 *    # COMPUTATION TIME: 105.915668 sec
 * 
 * 4. As output provide values of covariance matrix in form (put them to the report):
 *      cov( 1, 1)=      0.33330407
 *      cov( 2, 1)=      0.33330627
 *      cov( 2, 2)=        0.354141
 *      ....            ....
 *      cov(10, 9)=      0.30952735
 *      cov(10,10)=       0.14192182
 * 
 * 5. To speed up you work start with provided template below.
 * 
 * 6. Deadline for submission of report is: 09-03-2022
 * */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "mct_utils.h"

#define NELEMENTS 134217728

int main(int argc, char **argv)
{
    int i, j;

    // allocate memory
    double *var[5];
    for (i = 0; i < 5; i++)
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

    // generate additional random variables v_6, ..., v_10 and compute covariance matrix
    // make computation only for elements cov(i,j) where i<=j (see printing statment)
    double cov[10][10]; // storage for covaraince matrix
    b_t();              // start timing

    // TODO: add your code here!

    double tcmp = e_t(); // stop timing
    printf("# COMPUTATION TIME: %f sec\n", tcmp);

    // print results
    for (i = 0; i < 10; i++)
        for (j = 0; j <= i; j++)
            printf("cov(%2d,%2d)=%16.8g\n", i + 1, j + 1, cov[i][j]);

    return (EXIT_SUCCESS);
}
