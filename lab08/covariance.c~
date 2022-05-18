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
  FILE* fout = fopen("time_mpi.txt","a"); // file for storring timew
  if(fout == NULL) {
    printf("# Error: cannot open the file");
    return EXIT_FAILURE;
  }
  
  // Initialize MPI
  int ip, np; // basic MPI indicators
  char inmsg, outmsg='x';
  MPI_Status Stat; // required variable for receive routines
  MPI_Init( &argc, &argv ); // set up the parallel world
  MPI_Comm_size( MPI_COMM_WORLD, &np ); // total number of processes
  MPI_Comm_rank( MPI_COMM_WORLD, &ip ); // id of process st 0 <= ip < np

  // MPI I/O pointers
  MPI_File fh;
  MPI_Offset my_offset;
  MPI_Status status;
  int count;

  long int Nip, i0;  
  int i, j, in;
  double *var[10];
  /* for(i=0; i<5; i++) cppmallocl(var[i],NELEMENTS,double);      */
  FILE * fp;
  char file_name[128];
          

  // distribute elements over processe
  Nip = getNip(ip, np, NELEMENTS);
  i0 = get_i0(ip, np, NELEMENTS);
  size_t size = sizeof(double)*Nip;
  // printf("Process ip=%3d: Number of elements: %ld; Range: [%ld,%ld]; Size of array: %ld MB\n", ip, Nip, i0, i0+Nip-1, size/1024/1024);

  // Allocate memory
  for(i=0; i<10; i++) cppmallocl(var[i], Nip, double);
    
    
  // Read binary data
  MPI_Barrier(MPI_COMM_WORLD);
  // MPI_Rev(void* data, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm communicator, MPI_Status* status);
  if(ip>0){
    // printf("\t- process %3d waiting...", ip);
    MPI_Recv(&inmsg, 1, MPI_CHAR,ip-1, 1, MPI_COMM_WORLD, &Stat ); // wait for signal from previous process
    // printf("...received\n");
  }

  b_t(); // start timing
  printf("\t- process %3d loading data:\n", ip);
  for(i=0; i<5; i++){
      	
    sprintf(file_name, "/home2/archive/MCT-2022/lab1/var%d.dat", i+1);
    // printf("# READING DATA: `%s`...\n", file_name);
    MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    my_offset=(MPI_Offset)(i0*sizeof(double));
    MPI_File_seek(fh, my_offset, MPI_SEEK_SET);

    MPI_File_read_all(fh, var[i], Nip, MPI_DOUBLE, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &count);
    if(count!=Nip) { printf("# ERROR: Cannot read `%s`!\n", file_name); return EXIT_FAILURE; }

    // close file
    MPI_File_close(&fh);
	
  }// end of reading files for loop
  double read_time = e_t();

  b_t();
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
    
   
  // generate additional random variables v_6, ..., v_10 and compute covariance matrix
  // make computation only for elements cov(i,j) where i<=j (see printing statment)
  double cov[10][10] = { 0 }; // global storage for covaraince matrix
  double priv_cov[10][10] = { 0 }; // private storage for covariance for each process
  double avg_var[10] = { 0 }; // global average array
  double priv_avg[10] = { 0 }; // private average array (for each process)
  double var0i, var1i, var2i, var3i, var4i, var5i, var6i, var7i, var8i, var9i; //memory optimalization
  MPI_Barrier(MPI_COMM_WORLD);
  b_t();              // start timing

  //generate additional random variables
  for(in = 0; in < Nip; in++){ // calculating partial sum for a chunk
    /* priv_avg[0] += var[0][in]; */
    /* priv_avg[1] += var[1][in]; */
    /* priv_avg[2] += var[2][in]; */
    /* priv_avg[3] += var[3][in]; */
    /* priv_avg[4] += var[4][in]; */

    /* var[5][in] = sin(var[1][in]) + sin(var[0][in]); */
    /* var[6][in] = exp(var[2][in]) - exp(-1.*var[4][in]); */
    /* var[7][in] = sin(var[3][in])*cos(var[0][in]) + cos(var[3][in])*sin(var[2][in]); */
    /* var[8][in] = hypot(var[2][in], var[1][in]); */
    /* var[9][in] = cbrt(var[3][in]); */
    
    /* priv_avg[5] += var[5][in]; */
    /* priv_avg[6] += var[6][in]; */
    /* priv_avg[7] += var[7][in]; */
    /* priv_avg[8] += var[8][in]; */
    /* priv_avg[9] += var[9][in]; */
    
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
    // total of 5 + 15 = 20 (was 21+15 = 36)
  } 
  
  for(i = 0; i<10; i++) priv_avg[i] /=  NELEMENTS; // calculating partial average
  // MPI_Reduce(void* send_data, void* recv_data, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm communicatior);
  // MPI_Alleduce(void* send_data, void* recv_data, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm communicatior);
  // calculate gloabl average
  MPI_Allreduce(priv_avg, avg_var, 10, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

  // calculate covariance
  double tmp_cov = 0;
  for (i = 0; i < 10; i++){
    for(j = 0; j <= i; j++){
      tmp_cov = 0;
      for( in=0; in<Nip; in++)
	tmp_cov += (var[i][in] - avg_var[i]) * (var[j][in] - avg_var[j]);
      priv_cov[i][j] = tmp_cov / (NELEMENTS - 1);
    }
  }
  MPI_Reduce(priv_cov, cov, 100, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  double tcmp = e_t(); // stop timing

  // writing time
  b_t(); // start timing
  for(i=0; i<5; i++){      	
    sprintf(file_name, "/home/dteam201/MCT/lab08/var%d.dat", i+6);
    // printf("# READING DATA: `%s`...\n", file_name);
    MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    my_offset=(MPI_Offset)(i0*sizeof(double));
    MPI_File_seek(fh, my_offset, MPI_SEEK_SET);

    MPI_File_write_all(fh, var[i], Nip, MPI_DOUBLE, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &count);
    if(count!=Nip) { printf("# ERROR: Cannot write `%s`!\n", file_name); return EXIT_FAILURE; }

    // close file
    MPI_File_close(&fh);
	
  }// end of writing files for loop
  double write_time = e_t();
  printf("Reading time: &lf s", e_t);

  if(ip==0){
    printf("# READING TIME: %f sec\n", read_time);
    printf("# COMPUTATION TIME: %f sec\n", tcmp);
    printf("# WRITING TIME: %f sec\n", write_time);
    printf("# number of processes = %d\n",np);
    fprintf(fout,"%d\t%lf\t%lf\t%lf\n",np,read_time,tcmp,write_time);
    fclose(fout);
  }
  


  
  MPI_Barrier(MPI_COMM_WORLD);
  // print results
  if(ip==0){
    for(i=0; i<10; i++)
      for(j=0; j<=i; j++)
	printf("cov(%2d,%2d)=%16.8g\n", i+1, j+1, cov[i][j]);
  }

  // done with MPI
  MPI_Finalize();
      
  return( EXIT_SUCCESS ) ;
}

