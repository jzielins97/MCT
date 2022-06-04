/**
 * Template for lab9, "Modern Computing Technologies"
 * 
 * Author: M.Grunwald

* When compiling on dwarf:
 *  use computing node: ssh61 - ssh67
 *  module load cuda/9.0 
 * 
 * source scl_source enable devtoolset-7 python27 
 * nvcc -arch sm_35 -O3 laplace2d.cu -o laplace2d -lcudart -lcufft -lm

 * */
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/complex.h>

#include "mct_utils.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// See for reference:
// https://en.cppreference.com/w/c/numeric/complex
#include <complex.h>



//=========================================================
// Functions
//=========================================================
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
    for (ixyz = 0; ixyz < N; ixyz++)
    {
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

__global__ void calcFkxy(int nx, int ny, cufftDoubleComplex* fk, cufftDoubleComplex* fkxy)
{
    // Get our global thread ID
    size_t ixy = blockIdx.x*blockDim.x+threadIdx.x;
    cufftDoubleComplex z;
    
    double Lx = 100.0; // width in x-direction
    double Ly = 100.0; // width in y-direction

    double dx = Lx/nx;
    double dy = Ly/ny;

    int ix, iy;
    double kx, ky, k2, factor;

    // Make sure we do not go out of bounds
    if(ixy < nx*ny){
        ix = ixy/ny;
        iy = ixy -ix*ny;

        if(ix<nx/2) kx=2.*M_PI/(dx*nx)*(ix);
        else kx=2.*M_PI/(dx*nx)*(ix-nx);
                
        if(iy<ny/2) ky=2.*M_PI/(dy*ny)*(iy);
        else ky=2.*M_PI/(dy*ny)*(iy-ny);

        k2 = -1.0*(kx*kx + ky*ky);
        factor = k2/(nx*2*(ny - 1)); //because here ny is "in reality" ny/2 + 1, so we have to turn it back into "real" ny, so 2*(ny-1) 

        z = fk[ixy];
        z.x*=factor;
        z.y*=factor;
        fkxy[ixy]=z;
    }

}

//=========================================================
// Main
//=========================================================

int main()
{

    //----------------------------------
    // Settings
    //----------------------------------
    bool debug = false;
    int nx = 4112; // number of points in x-direction
    int ny = 4112; // number of points in y-direction
    int blockSize=1024; // Number of threads in each thread block
    int gridSize; // Number of thread blocks in grid

    double Lx = 100.0; // width in x-direction
    double Ly = 100.0; // width in y-direction

    double x0 = -Lx / 2;
    double y0 = -Ly / 2;

    double dx = Lx/nx;
    double dy = Ly/ny;

    size_t sizeReal =  nx*ny*sizeof(double);
    size_t sizeComplex =  nx*(ny/2 +1)*sizeof(std::complex<double>);
 

    double *fxy_host; // function
    double *fxy_device;  
    std::complex<double> *fk_device; // transformed function
    std::complex<double> *fkxy_device;
    double *laplacefxy_host; // computed numerically
    double *laplacefxy_device;
    double *formula_laplacefxy; // computed according formula

    cudaSetDevice(0); // to see devices use: nvidia-smi


    cppmallocl(fxy_host, nx*ny, double);
    cudaMalloc(&fxy_device, sizeReal);
    cudaMalloc(&fk_device, sizeComplex);
    cudaMalloc(&fkxy_device, sizeComplex);
    cppmallocl(laplacefxy_host, nx*ny, double);
    cudaMalloc(&laplacefxy_device, sizeReal);
    cppmallocl(formula_laplacefxy, nx*ny, double);

    //------------------------------------------------
    // Initialising function(s), copying to the device
    //------------------------------------------------

    int ix, iy, ixy=0;

    b_t(); //start timing (for initializing)

    for (ix = 0; ix < nx; ix++){
    	for (iy = 0; iy < ny; iy++){
	    // function
            fxy_host[ixy] = function_xy(x0 + dx * ix, y0 + dy * iy);

            // result for comparion
            formula_laplacefxy[ixy] = laplace_function_xy(x0 + dx * ix, y0 + dy * iy);
	    ixy++;
	}
    }
    double tInit = e_t(); // stop timing
    printf("# FUNCTIONS INITIALIZING TIME: %f sec\n", tInit);

    b_t(); //start timing (for sending)
    cudaMemcpy(fxy_device, fxy_host, sizeReal, cudaMemcpyHostToDevice);
    double tSend = e_t(); // stop timing
    printf("# COPYING (HOST->DEVICE) TIME: %f sec\n", tSend);

	//-------------------------------------------------------------------
	// preparing plans, transforming forward
	//--------------------------------------------------------------------
    cufftHandle plan_f, plan_b;
    cufftResult result;
    size_t required_size_f, required_size_b;

    b_t(); //start timing (for computing)
    result = cufftPlan2d(&plan_f, nx, ny, CUFFT_D2Z);
    if(debug)printf("Handle forward allocation - %d\n",(int) result);

    result = cufftEstimate2d(nx, ny, CUFFT_D2Z, &required_size_f);
    if(debug)printf("Worksize forward allocation - %d, size: %zd\n",(int) result, required_size_f);

    result = cufftMakePlan2d(plan_f, nx, ny, CUFFT_D2Z, &required_size_f);
    if(debug)printf("Plan forward making - %d\n",(int) result);

    result = cufftExecD2Z(plan_f, (cufftDoubleReal*)fxy_device, (cufftDoubleComplex*)fk_device);
    if(debug)printf("Plan forward executing - %d\n",(int) result);

	//----------------------------------------------------------
	// creating transformed laplace
	//---------------------------------------------------------

    gridSize = (int)ceil((float)nx*(ny/2+1)/blockSize);
    calcFkxy<<<gridSize, blockSize>>>(nx, ny/2+1, (cufftDoubleComplex*)fk_device, (cufftDoubleComplex*)fkxy_device);

	//--------------------------------------------------
	// transforming backwards, coping to the host
	//--------------------------------------------------
    result = cufftPlan2d(&plan_b, nx, ny, CUFFT_Z2D);
    if(debug)printf("Handle backwards allocation - %d\n",(int) result);

    result = cufftEstimate2d(nx, ny, CUFFT_Z2D, &required_size_b);
    if(debug)printf("Worksize backwards allocation - %d, size: %zd\n",(int) result, required_size_b);

    result = cufftMakePlan2d(plan_b, nx, ny, CUFFT_Z2D, &required_size_b);
    if(debug)printf("Plan backwards making - %d\n",(int) result);

    result = cufftExecZ2D(plan_b, (cufftDoubleComplex*)fkxy_device, (cufftDoubleReal*)laplacefxy_device);
    if(debug)printf("Plan backwards executing - %d\n",(int) result);

    double tCmpt = e_t(); // stop timing
    printf("# COMPUTING (FORWARD->KERNEL->BACKWARD) TIME: %f sec\n", tCmpt);

    b_t();//start timing (receive)
    cudaMemcpy(laplacefxy_host, laplacefxy_device, sizeReal, cudaMemcpyDeviceToHost);
    double tRecv = e_t(); // stop timing
    printf("# COPYING (DEVICE->HOST) TIME: %f sec\n", tRecv);
    printf("# TOTAL TIME: %f sec\n", tInit+tSend+tCmpt+tRecv);

	//--------------------------------------------------
	// checking results
	//--------------------------------------------------
    test_array_diff(nx * ny, laplacefxy_host, formula_laplacefxy);

    if(debug){
        ixy=0;
        for(ix=0; ix<nx; ix++){ 
            for(iy=0; iy<ny; iy++){    
                if(ix==nx/2) printf("%10.6g %10.6f %10.6f %10.6f\n", x0 + dx*ix, y0 + dy*iy, laplacefxy_host[ixy], formula_laplacefxy[ixy]);
                ixy++;
            }
        }        
    }


    cudaFree(fxy_device);
    cudaFree(fk_device);
    cudaFree(fkxy_device);
    cudaFree(laplacefxy_device);

    free(fxy_host);
    free(laplacefxy_host);
    free(formula_laplacefxy);

    return 1;
}
