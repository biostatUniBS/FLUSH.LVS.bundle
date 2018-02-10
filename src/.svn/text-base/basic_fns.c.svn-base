/*********************************************************************
 **
 ** file: basic_fns.c
 **
 ** Aim: Function used by lvs_rlm.c
 **
 ** Copyright (C) 2007 Stefano Calza 
 ** Mainly adapted from code by Ben M .Bolstad <bolstad@stat.berkeley.edu> 2003
 **
 ** 
 ** All code was written by Ben M .Bolstad <bolstad@stat.berkeley.edu> 2003-2007
 ** 
 ** created on: May 10, 2007
 **
 ** Last modified: May 10, 2007
 **
 ** The aim is to provide a set of functions to fit a 
 ** robust linear models to affy data in order to get arrays & probes variance components. 
 ** It is focused on huber regression. Code is inspired by rlm() method which is
 ** part of the MASS package bundle and from affyPLM set of functions.
 */

#include "basic_fns.h"

/* ok */
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



/**********************************************************
 **
 ** int sort_double(const void *a1,const void *a2)
 ** 
 ** a comparison function used when sorting doubles.
 **
 **********************************************************/

/*
int lvs_sort_double(const double *a1,const double *a2){
  if (*a1 < *a2)
    return (-1);
  if (*a1 > *a2)
    return (1);
  return 0;
}
*/


/**************************************************************************
 **
 ** double median(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the median of *x
 **
 *************************************************************************/

double  lvs_median(double *x, int length){
  int i;
  int half;
  double med;
  double *buffer = Calloc(length,double);
  
  memcpy(buffer,x,length*sizeof(double));

  half = (length + 1)/2;

  rPsort(buffer, length, half-1);
  med = buffer[half-1];
  if (length % 2 == 0){
    rPsort(buffer, length, half);
    med = (med + buffer[half])/2.0;
  }
  
  Free(buffer);
  return med;
}


/**************************************************************************
 **
 ** double median_nocopy(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the median of *x. note x is not order preserved when this function
 ** is called.
 **
 *************************************************************************/
/*
double  lvs_median_nocopy(double *x, int length){
  int i;
  int half;
  double med;
  double *buffer = x;
  
  memcpy(buffer,x,length*sizeof(double));
  
  half = (length + 1)/2;
  
  rPsort(buffer, length, half-1);
  med = buffer[half-1];
  if (length % 2 == 0){
    rPsort(buffer, length, half);
    med = (med + buffer[half])/2.0;
  }
  
 
  return med;
}
*/


/***************************************************************
 **
 ** double irls_delta(double *old, double *new, int length)
 **
 ** double *old - previous value of a vector
 ** double *new - new value of a vector
 ** int length - length of vector
 **
 ** this function computes the sum of the difference of two vectors
 ** divides this by the sum squared of the old datavector.
 **
 ** the aim of this function is compute something to test for 
 ** convergence in the iteratively reweighted least squares (IRLS)
 ** 
 **
 ** Author: Ben Bolstad
 **
 **************************************************************/

double lvs_irls_delta(double *old, double *new, int length){
  int i=0;
  double sum = 0.0;
  double sum2 =0.0;
  double divisor=1e-20;

  for (i=0; i < length; i++){
    sum = sum + (old[i] - new[i])*(old[i]-new[i]);
    sum2 = sum2 + old[i]*old[i];
  }
  
  if(sum2 >= divisor){
    divisor = sum2;
  }

  return sqrt(sum/divisor); 
} 


/**********************************************************************************
 **
 ** double med_abs(double *x, int length)
 **
 ** double *x - a vector of data
 ** int length - length of the vector.
 ** 
 ** returns the median of the absolute values.
 **
 ** computes the median of the absolute values of a given vector.
 **
 ** Author: Ben Bolstad
 **********************************************************************************/

double lvs_med_abs(double *x, int length){
  int i;
  double med_abs;
  double *buffer = Calloc(length,double);

  for (i = 0; i < length; i++)
    buffer[i] = fabs(x[i]);
  
  med_abs = lvs_median(buffer,length);
    
  Free(buffer);
  return(med_abs);
}


double lvs_psi_huber(double u, double k,int deriv)
{
  
  if (deriv == 0){
    if ( 1 < k/fabs(u)){
      return 1.0;
    } else {
      return  k/fabs(u);
    }
  } else if (deriv == 1){
    if (fabs(u) <= k){
      return 1.0;
    } else {
      return 0.0;
    }
  } else {
    if (fabs(u) <= k){
      return u;
    } else {
      if (u < 0){
	return -k;
      } else {
	return k;
      }
    }
  }
}
