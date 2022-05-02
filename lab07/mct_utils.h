/**
 * Modern Computing Technologies, WUT, 2021
 * 
 * This file contains useful snippets of a code
 * */ 

/// allocation of memory, local verions, not involving MPI
#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "error: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    } 
    

#include <sys/time.h>
#include <time.h> /* for ctime() */
static double t_gettimeofday ;
static struct timeval s ;

void b_t( void ) 
{ /* hack together a clock w/ microsecond resolution */
  gettimeofday( &s , NULL ) ;
  t_gettimeofday = s.tv_sec + 1e-6 * s.tv_usec ;
}

double e_t() 
{
    gettimeofday( &s , NULL ) ;
    return s.tv_sec + 1e-6 * s.tv_usec - t_gettimeofday ;
}


/**
 * Function for uniform distribution of array elemnts
 * @param ip index of MPI process
 * @param np number of MPI processes 
 * @param N number of array elements
 * @return chunk size
 * */
long int getNip( int ip , int np , long int N ) 
{
  long int Nip;
  Nip  = N / np ;
  if ( ip < (N % np) ) Nip++ ;
  return Nip;
}

/**
 * Function for uniform distribution of array elemnts
 * @param ip index of MPI process
 * @param np number of MPI processes 
 * @param N number of array elements
 * @return index of first element in global array. 
 * */
long int get_i0( int ip , int np , long int N ) 
{
  long int Nip;
  long int i0=0; 
  int i;
  for(i=0; i<ip; i++)
    {
      Nip  = N / np ;
      if ( i < (N % np) ) Nip++ ;
      i0+=Nip;
    }
  return i0;
}

