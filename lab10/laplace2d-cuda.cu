
/**
 * Task8, "Modern Computing Technologies"
 * 
 * Author: Jakub Zielinski
 * Date: 5.06.2022
 * 
 * nvcc -arch sm_35 -O3 laplace2d-cuda.cu -o laplace2d-cuda -lcudart -lcufft -lm
 * */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "../mct_utils.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define BATCH 1

// See for reference:
// https://en.cppreference.com/w/c/numeric/complex
#include <complex.h>

// See for reference:
// https://docs.nvidia.com/cuda/cufft/index.html
#include <cufft.h>

/**
 * Test function
 * */
double function_xy(double x, double y)
{
#define A -0.03
#define B -0.01
#define C -0.005
  return exp(A * x * x + B * y * y + C * x * y);
}

/**
 * Analytical result
 * Use it for checking correctness!
 * */
double laplace_function_xy(double x, double y)
{
  double rdf2d_dx = function_xy(x, y) * (2. * A * x + C * y); // d/dx
  double rdf2d_dy = function_xy(x, y) * (2. * B * y + C * x); // d/dy
  
  double rlaplacef2d = rdf2d_dx * (2. * A * x + C * y) + function_xy(x, y) * (2. * A) + rdf2d_dy * (2. * B * y + C * x) + function_xy(x, y) * (2. * B); // laplace
  
  return rlaplacef2d;

#undef A
#undef B
#undef C
}

/**
 * You can use this function to check diff between two arrays
 * */
void test_array_diff(int N, double *a, double *b)
{
  int ixyz = 0;
  
  double d, d2;
  double maxd2 = 0.0;
  double sumd2 = 0.0;
  for (ixyz = 0; ixyz < N; ixyz++){
    d = a[ixyz] - b[ixyz];
    
    d2 = d * d;
    sumd2 += d2;
    if (d2 > maxd2)
      maxd2 = d2;
  }

  printf("#    COMPARISON RESULTS:\n");
  printf("#           |max[a-b]| : %16.8g\n", sqrt(maxd2));
  printf("#         SUM[(a-b)^2] : %16.8g\n", sumd2);
  printf("# SQRT(SUM[(a-b)^2])/N : %16.8g\n", sqrt(sumd2) / N);
}


// CUDA kernels -------------------------->
__global__ void calculate_Fk(cufftDoubleComplex* fK, int nx, int ny, double dx, double dy)
{
  cufftDoubleComplex z;
  size_t ixy = blockIdx.x*blockDim.x+threadIdx.x;
  int ny_prim = ny/2+1;
  int ix = ixy/ny_prim;
  int iy = ixy - ix*ny_prim;
  double kx = 0;
  double ky = 0;
  if(ixy < nx * (ny/2+1) ){  
    if(ix<nx/2) kx = 2.*M_PI/(dx*nx)*(ix   );
    else        kx = 2.*M_PI/(dx*nx)*(ix-nx);
      
    if(iy<ny/2) ky = 2.*M_PI/(dy*ny)*(iy   );
    else        ky = 2.*M_PI/(dy*ny)*(iy-ny);

  //   // recalculate fK real and imaginary part, as cufftComplex  doesn't support some calculations 
    z=fK[ixy];
    z.x *= (-kx*kx - ky*ky) / (nx*ny);
    z.y *= (-kx*kx - ky*ky) / (nx*ny);
    fK[ixy] = z;
  }
}
//----------------------------------------<

int main()
{
  double GB = 1. /1024/1024/1024; // conversion to GB
  // plans for the FFTW
  cufftHandle plan_f;
  cufftHandle plan_b;
  
  // Settings
  int nx = 4112; // number of points in x-direction
  int ny = 4112; // number of points in y-direction	
  // int n_r[2] = {nx,ny}; // dimensions for the fxy
  // int n_c[2] = {nx,ny/2+1}; // dimensions for the fK
  
  double Lx = 100.0; // width in x-direction
  double Ly = 100.0; // width in y-direction
  
  double x0 = -Lx / 2;
  double y0 = -Ly / 2;
  
  double dx = Lx / nx;
  double dy = Ly / ny;
  
  double *h_fxy;
  double *formula_laplacefxy;

  cufftDoubleReal *fxy;                // function
  cufftDoubleComplex *fK;              // array with f(kx, ky)

  double init_t, send_t, calc_t, recv_t; // for timing

  int ix, iy, ixy; // for iterating
  int blockSize, gridSize;

  size_t realSize = nx * ny * sizeof(cufftDoubleReal);
  size_t complexSize = nx * (ny/2+1) * sizeof(cufftDoubleComplex);

  cudaSetDevice(0); // select GPU for CUDA
  
  printf("Allocate memory\n");
  //---- Allocate memory on the host ------>
  // printf("\tAllocate memory for fxy at host\n");
  cppmallocl(h_fxy, nx * ny, double);
  // printf("\tAllocate memory for fromula_laplacefxy at host\n");
  cppmallocl(formula_laplacefxy, nx * ny, double); // this doesn't have to be coppied to CUDA
  //---------------------------------------<
  
  //---- Allocate memory on the device ---->
  // printf("\tAllocate memory for fxy at device (cuda)\n");
  cudaMalloc(&fxy, realSize);
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to allocate\n");
    return -1;
  }
  // printf("\tAllocate memory for fK at device (cuda)\n");
  cudaMalloc(&fK, complexSize); // this doesn't have to saved on the host
  if (cudaGetLastError() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to allocate\n");
    return -1;
  }  
  //---------------------------------------<


  
  //--- initialize vectors on the host ---->
  printf("Initializing the functions vectors\n");
  b_t();
  ixy = 0;
  for (ix = 0; ix < nx; ix++){
    for (iy = 0; iy < ny; iy++){
      // function
      h_fxy[ixy] = function_xy(x0 + dx * ix, y0 + dy * iy);
      
      // result for comparion
      formula_laplacefxy[ixy] = laplace_function_xy(x0 + dx * ix, y0 + dy * iy);
      ixy++;
    }
  }
  init_t = e_t();
  printf("Initialization time of the fxy vector: %lf s\n", init_t);
  //---------------------------------------<
  
  //---- Copy host vectors to device ------>
  b_t();
  cudaMemcpy( fxy, h_fxy, realSize, cudaMemcpyHostToDevice);
  send_t = e_t();
  printf("Copy to device time: %f [sec]; bandwidth=%f [GB/sec]\n", send_t, sizeof(double)*nx*ny*GB/send_t);
  //---------------------------------------<

  /********* start the calculation ********/
  b_t();
  //------- create the forward plan ------->
  if(cufftPlan2d(&plan_f, nx, ny, CUFFT_D2Z) != CUFFT_SUCCESS){
    fprintf(stderr,"CUFFT error: unable to create forward plan\n");
  }
  //--------------------------------------<
  
  //------- execute the forward plan ------>
  if(cufftExecD2Z(plan_f, fxy, fK) != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: unable to execute the forward transform\n");
  }//--------------------------------------<

  //------------ construct F(kx,ky) ------->
  blockSize = 1024; // Number of threads in each thread block
  gridSize = (int)ceil((float)nx*(ny/2+1)/blockSize); // Number of thread blocks in grid
 
  // Execute the kernel
  calculate_Fk<<<gridSize, blockSize>>>(fK, nx, ny, dx, dy);
  // calcFkxy<<<gridSize, blockSize>>>(nx, ny/2+1, (cufftDoubleComplex*)fK);
  //---------------------------------------<

  //------ create the backward plan ------->
  if(cufftPlan2d(&plan_b, nx, ny, CUFFT_Z2D) != CUFFT_SUCCESS){
    fprintf(stderr,"CUFFT error: unable to create handle for backward plan\n");
  }
  //--------------------------------------<
  
  //----- execute the backward plan ------->
  if(cufftExecZ2D(plan_b, fK, fxy) != CUFFT_SUCCESS){
    fprintf(stderr, "CUFFT error: unable to execute the backward transform\n");
  }//--------------------------------------<
  /********* end of the calculation *******/
  calc_t = e_t();
  printf("Computation time: %lf s\n", calc_t);
    
  //------- Copy results to host ---------->
  b_t();
  cudaMemcpy( h_fxy, fxy, realSize, cudaMemcpyDeviceToHost);
  recv_t = e_t();
  printf("Copy to host time: %f [sec]; bandwidth=%f [GB/sec]\n", recv_t, sizeof(double)*nx*ny*GB/recv_t);
  //---------------------------------------<

  // Check correctness of computation
  test_array_diff(nx * ny, h_fxy, formula_laplacefxy);
  
  // FILE* fout = fopen("result.dat","w");
  // ixy=0;
  // for(ix=0; ix<nx; ix++){
  //   for(iy=0; iy<ny; iy++)
  // 	{
  // 	  if(ix==nx/2) fprintf(fout,"%10.6f %10.6f %10.6f %10.6f\n", x0 + dx*ix, y0 + dy*iy, h_fxy[ixy], formula_laplacefxy[ixy]);
  
  // 	  ixy++;
  // 	}
  // }
  // fclose(fout);


  cufftDestroy(plan_f);
  cufftDestroy(plan_b);
  cudaFree(fxy);
  cudaFree(fK);
  

  return 1;
}
