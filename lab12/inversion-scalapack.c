/**
 * Author: Gabriel Wlazłowski
 * Example code for "Modern Computing Technologies"
 * 
 * mpicc inversion-scalapack.c -o inversion-scalapack -L/usr/lib64/openmpi/lib -lscalapack -L/usr/lib64/atlas/ -llapack 
 * */ 


#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>

#include <mpi.h>

#include "mct_utils.h"

// BLACS
extern void blacs_get_(int* ICONTXT, int *WHAT, int *VAL);
extern void blacs_gridinit_(int* ICONTXT, char* ORDER, int* NPROW, int* NPCOL );
extern void blacs_gridinfo_(int* ICONTXT, int* NPROW, int* NPCOL, int* MYPROW, int* MYPCOL );

// BLOCK-CYCLIC DISTRIBUTION MATRIX MAPPING FUNCTIONS
// get number of rows/columns of submatrix
// http://www.netlib.org/scalapack/explore-html/d4/d48/numroc_8f_source.html
extern int numroc_(int * M, int * MB, int * ROWID, int * ISRCPROC, int * NPROW );

// convert global index gi into local li (global to local)
// http://www.netlib.org/scalapack/explore-html/d1/daa/indxg2l_8f_source.html
extern int indxg2l_(int * INDXGLOB, int * NB, int * IPROC, int * ISRCPROC, int * NPROCS );

// convert global index gi into process coordinate (global to process)
// http://www.netlib.org/scalapack/explore-html/d6/d88/indxg2p_8f_source.html
extern int indxg2p_(int * INDXGLOB, int * NB, int * IPROC, int * ISRCPROC, int * NPROCS );

// convert local index li into global gi (local to global)
// http://www.netlib.org/scalapack/explore-html/d4/deb/indxl2g_8f_source.html
extern int indxl2g_(int * INDXLOC, int * NB, int * IPROC, int * ISRCPROC, int * NPROCS );

// function for creating matrix description
// http://www.netlib.org/scalapack/explore-html/dd/d22/descinit_8f_source.html
extern void descinit_(int * DESC, int * M, int * N, int * MB, int * NB, 
                      int * IRSRC, int * ICSRC, int * ICTXT, int * LLD, int *INFO );

// ScaLAPACK
// http://www.netlib.org/scalapack/explore-html/df/dfe/pdgetrf_8f_source.html
extern void pdgetrf_(int * M, int * N, double * A, int * IA, int * JA, int * DESCA, int * IPIV, int * INFO );

// http://www.netlib.org/scalapack/explore-html/d3/df3/pdgetri_8f.html
extern void pdgetri_(int * N, double * A, int * IA, int * JA, int * DESCA, int * IPIV, double * WORK, int * LWORK, int * IWORK, int * LIWORK, int * INFO);

// http://www.netlib.org/scalapack/explore-html/d6/da2/pdgemm___8c_source.html
extern void pdgemm_(char * TRANSA, char * TRANSB,
               int * M, int * N, int * K,
               double * ALPHA,
               double * A, int * IA, int * JA, int * DESCA,
               double * B, int * IB, int * JB, int * DESCB,
               double * BETA,
               double * C, int * IC, int * JC, int * DESCC );

// matrix size [MxN=M]
#define M 4000
// size of block in M direction 
#define MB 8
// size of block in N direction
#define NB 8
/**
 * Function returns value of matrix element A_{ij}
 * i,j iterates from 1 to M as in standard mathematical notation
 * */
double matrix_element(int i, int j)
{
    // random matrix element
    return (double)rand() / (double)RAND_MAX ;
}

// translate matrix element into index into array, column major 
#define IDX(i,j,n) (j-1)*n + (i-1)

int main(int argc, char *argv[]) 
{
    
    // Initilize MPI 
    int ip, np; // basic MPI indicators
    MPI_Init( &argc , &argv ) ; /* set up the parallel WORLD */
    MPI_Comm_size( MPI_COMM_WORLD , &np ) ; /* total number of processes */
    MPI_Comm_rank( MPI_COMM_WORLD , &ip ) ; /* id of process st 0 <= ip < np */    
    
    srand(ip+1); // initilize random number generator diffrently for each process
    
    if(ip==0)
        printf("MPIRUN EXECUTED WITH np=%d PROCESSES\n", np);
        
    int dims[2] = {0,0};
    MPI_Dims_create(np, 2, dims); // however, you can also set nprow and npcol by hand, keeping constraing nprow*npcol=np
    int nprow = dims[0]; // cartesian direction 0
    int npcol = dims[1]; // cartesian direction 1
    
    if(ip==0)
        printf("SETTING BLACS GRID OF SIZE [%d x %d]\n", nprow, npcol);

    // Create BLACS grid
    // Get a default BLACS context
    int zero=0, mone=-1, one=1;
    int context;
    blacs_get_(&mone, &zero, &context);
    
    // Initialize the BLACS context
    int iprow, ipcol;
    blacs_gridinit_(&context, "R", &nprow, &npcol);
    blacs_gridinfo_(&context, &nprow, &npcol, &iprow, &ipcol);
    
    if(iprow==-1 || ipcol==-1)
    {   
        printf("ERROR: ip=%d CANNOT CREATE BLACS GRID!\n", ip);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if("BLACS GRID CREATED\n");
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Determine size of local storage for each process
    double *A, *B; // matrix
    int gM = M; // number of rows in global matrix
    int gN = M; // number of columns in global matrix
    
    int mb = MB; // block size - parameter
    int nb = NB; // block size - parameter
    int lM = numroc_(&gM, &mb, &iprow, &zero, &nprow ); // number of rows for local matrix
    int lN = numroc_(&gN, &nb, &ipcol, &zero, &npcol ); // number of columns for local matrix
    
    // and allocate memory
    printf("PROCESS %6d -> [%3d,%3d] ALLOCATES SUBMATRIX OF SIZE [%d x %d]\n", ip, ipcol, iprow, lM, lN);
    cppmallocl(A,lM*lN,double);
    cppmallocl(B,lM*lN,double); // storage for copy 
    
    // fill matrix with values
    if(ip==0) printf("INITIALIZE MATRIX...\n");
    int li, lj; // local indices
    int gi, gj;
    for(lj=1; lj<=lN; lj++) // for each (sub)column  < -- note column major !!!
        for(li=1; li<=lM; li++) // for each (sub)row
        {
            // map local index (li,lj) into global one (gi,gj)
            gi = indxl2g_(&li, &mb, &iprow, &zero, &nprow ); // for row
            gj = indxl2g_(&lj, &nb, &ipcol, &zero, &npcol ); // for column
            A[IDX(li,lj,lM)] = matrix_element(gi, gj);
        }
        
    // make copy of A 
    for(li=0; li<lM*lN; li++) B[li]=A[li];
        
    // Create matrix descriptor - required by PBLAS and ScaLAPACK routines
    if(ip==0) printf("CREATE MATRIX DESCRIPTOR...\n");
    int descA[ 9 ]; // decriptor object
    int info;
    descinit_(descA, &gM, &gN, &mb, &nb, &zero, &zero, &context, &lM, &info);    
    if(info!=0)   
    {
        printf("ERROR: ip=%d CANNOT MATRIX DECRIPTOR!\n", ip);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // invert matrix
    
    // memory query 
    // * IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
    if(ip==0) printf("ALLOCATING WORKSPACE...\n");
    int *ipiv;
    cppmallocl(ipiv,(gM+mb),int);
    int lwork=-1, liwork=-1;
    double testwork[1];    
    int testiwork[1];
    pdgetri_(&gN, A, &one, &one, descA, ipiv, testwork, &lwork, testiwork, &liwork, &info);
    if(info!=0)   
    {
        printf("ERROR: ip=%d CANNOT QUERY pdgetri! INFO=%d\n", ip, info);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    lwork=(int)testwork[0];
    double *work;
    cppmallocl(work,lwork,double);
    liwork=(int)testiwork[0];
    int *iwork;
    cppmallocl(iwork,liwork,int);
    

    // do inversion
    MPI_Barrier(MPI_COMM_WORLD);
    b_t();
    
    if(ip==0) printf("INVERTING MATRIX...\n");
    pdgetrf_(&gM, &gN, A, &one, &one, descA, ipiv, &info );
    if(info!=0)   
    {
        printf("ERROR: ip=%d CANNOT EXECUTE pdgetrf! INFO=%d\n", ip, info);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    pdgetri_(&gN, A, &one, &one, descA, ipiv, work, &lwork, iwork, &liwork, &info);
    if(info!=0)   
    {
        printf("ERROR: ip=%d CANNOT EXECUTE pdgetri! INFO=%d\n", ip, info);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    double rt = e_t();
    if(ip==0) printf("INVERTING MATRIX TOOK %.2f SEC\n", rt);
    
    
    // Check correctness of inversion
    MPI_Barrier(MPI_COMM_WORLD);
    b_t();  
    
    // C = A * B shoule be equal ONE, if A is indeed inverted matrix
    double *C; // storage for result
    cppmallocl(C,lM*lN,double);
    double done=1.0, dzero=0.0;
    
    pdgemm_("N", "N", &gM, &gN, &gM, 
            &done, // alpha
            A, &one, &one, descA, // matrix A
            B, &one, &one, descA, // matrix B
            &dzero, // beta
            C, &one, &one, descA // matrix C = A*B
            );

    rt = e_t();
    if(ip==0) printf("MATRIX MULTIPLICATION TOOK %.2f SEC\n", rt);
    
    // check correctness 
    double epsilon=1.0e-9;
    for(lj=1; lj<=lN; lj++) // for each (sub)column  < -- note column major !!!
        for(li=1; li<=lM; li++) // for each (sub)row
        {
            // map local index (li,lj) into global one (gi,gj)
            gi = indxl2g_(&li, &mb, &iprow, &zero, &nprow ); // for row
            gj = indxl2g_(&lj, &nb, &ipcol, &zero, &npcol ); // for column
            
            if(gi==gj && fabs(C[IDX(li,lj,lM)]-1.0)>epsilon) 
                printf("PROBLEM WITH C[%d,%d]=%f\n", gi, gj, C[IDX(li,lj,lM)]);
            if(gi!=gj && fabs(C[IDX(li,lj,lM)]-0.0)>epsilon) 
                printf("PROBLEM WITH C[%d,%d]=%f\n", gi, gj, C[IDX(li,lj,lM)]);
        }
       
    MPI_Barrier(MPI_COMM_WORLD);
    if(ip==0) printf("TEST DONE!\n");
        
    // done with MPI  
    MPI_Finalize();
    
    return 0;
}
