/**
 * Modern Computing Technologies, WUT, 2022
 * 
 * Author: Jakub Zielinski
 * Date: 16.03.2022
 * 
 * Compilation command:
 *      mpicc covariance.c -o cov_mpi -lm 
 * Run command:
 *      mpirun -np 20 ./cov_mpi
 * */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

//include MPI library
#include "mpi.h" 

#include "mct_utils.h"

#define NELEMENTS 134217728

int main( int argc , char ** argv ) 
{
  // Initialize MPI
  int ip, np; // basic MPI indicators
  char inmsg, outmsg='x';
  MPI_Status Stat; // required variable for receive routines
  MPI_Init( &argc, &argv ); // set up the parallel world
  MPI_Comm_size( MPI_COMM_WORLD, &np ); // total number of processes
  MPI_Comm_rank( MPI_COMM_WORLD, &ip ); // id of process st 0 <= ip < np

  long int Nip, i0;  
  int i, j;
  double *var[5];
  /* for(i=0; i<5; i++) cppmallocl(var[i],NELEMENTS,double);      */
  FILE * fp;
  char file_name[128];
          

  // distribute elements over processe
  Nip = getNip(ip, np, NELEMENTS);
  i0 = get_i0(ip, np, NELEMENTS);
  size_t size = sizeof(double)*Nip;
  printf("Process ip=%3d: Number of elements: %ld; Range: [%ld,%ld]; Size of array: %ld MB\n", ip, Nip, i0, i0+Nip-1, size/1024/1024);

  // Allocate memory
  for(i=0; i<5; i++) cppmallocl(var[i], Nip, double);
    
    
  // Read binary data
  MPI_Barrier(MPI_COMM_WORLD);
  // MPI_Rev(void* data, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm communicator, MPI_Status* status);
  if(ip>0){
    // printf("\t- process %3d waiting...", ip);
    MPI_Recv(&inmsg, 1, MPI_CHAR,ip-1, 1, MPI_COMM_WORLD, &Stat ); // wait for signal from previous process
    // printf("...received\n");
  }

  b_t(); // start timing
  // printf("\t- process %3d loading data:\n", ip);
  for(i=0; i<5; i++){
      	
    sprintf(file_name, "/home2/archive/MCT-2022/lab1/var%d.dat", i+1);
    // printf("# READING DATA: `%s`...\n", file_name);
      
    fp = fopen(file_name, "rb");
    if(fp==NULL) { printf("# ERROR: Cannot open `%s`!\n", file_name); return EXIT_FAILURE; }
    
    //  move file pointer to the beginning of the chunk for the given process
    fseek( fp, i0*sizeof(double), SEEK_SET);
    // read chunk from the file
    size_t ftest = fread(var[i], sizeof(double), Nip, fp);
    // printf("\tRead from file #%1d: %ld; Range: [%ld,%ld]; Size of array: %ld MB\n", i, ftest, i0, i0+Nip-1, sizeof(double)*ftest/1024/1024);
    if(ftest!=Nip) { printf("# ERROR: Cannot read `%s`!\n", file_name); return EXIT_FAILURE; }

    // close file
    fclose(fp);
	
  }// end of reading files for loop

  if(ip!=(np-1)){
    //MPI_Send(void* data, int count, MPI_Datatype datatype, int destination, int tag, MPI_Comm communicator);
    // printf("\t- process %3d sending message to %3d...",ip, ip+1);
    MPI_Send(&outmsg, 1, MPI_CHAR, ip+1, 1, MPI_COMM_WORLD); // send signal to process ip+1 that reading of my chunks is finished
    // printf("sent\n");
  }// end of if(ip!=(np-1))
    
  double tio = e_t(); // stop timing
  if(ip==0) printf("# READ TIME: %f sec\n", tio);
    
  MPI_Barrier(MPI_COMM_WORLD);

  long int idx[10] = {10,13421673,25501328,41606496,53677091,73818750,83214911,93952210,106032910,132875451};
  for(i=0; i<5; i++) 
    for(j=0; j<10; j++)
      if(i0 <= idx[j] && i0+Nip > idx[j]){ printf("var[%d][%10ld]=%16.8g\n", i,idx[j],var[i][idx[j]-i0]); }
    
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


  // done with MPI
  MPI_Finalize();
      
  return( EXIT_SUCCESS ) ;
}

