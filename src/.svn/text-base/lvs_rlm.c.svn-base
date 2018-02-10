/*********************************************************************
 **
 ** file: lvs_rlm.c
 **
 ** Aim: implement robust linear models to be used for LVS normalization
 **
 ** Copyright (C) 2007 Stefano Calza 
 ** Mainly adapted from code by Ben M .Bolstad <bolstad@stat.berkeley.edu> 2003
 **
 ** created by Ben M .Bolstad <bolstad@stat.berkeley.edu> 2003-2007
 ** modified by: Stefano Calza <calza@med.unibs.it> 2007
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
#include "externals.h"

#include <R.h>
#include <Rmath.h>




void lvs_lm_wfit(double *x, double *y, double *w, int rows, int cols, double tol, 
		 double *out_beta, double *out_resids){
  int i,j;
  int ny = 1;
  int k;
  int numzero_weights = 0,nrows,totnumzero=0;
  
  double fittedvalue;

  double *wts = Calloc(rows,double);
  double *x_wts_f = Calloc(rows*cols,double);
  double *y_wts_f = Calloc(rows,double);
  double *beta = Calloc(cols,double);
  double *resid = Calloc(rows,double);
  double *qraux = Calloc(cols,double);
  double *qty = Calloc(rows,double);
  double *work = Calloc(2*cols,double);
  
  int *jpvt = Calloc(cols,int);

  for (i=0; i < rows; i++){
    if (w[i] == 0.0){
      totnumzero++;
    }
  }

  if (totnumzero > 0){
    /* we need to chop up the X and Y matricies a little more then removing the 
       observations that have weight = 0. In particular fit the model removing the weight
       0 observations. then to compute the residuals for the weight 0 observations used fitted
       values and observed X to compute fitted values */
    numzero_weights = 0;
    for (i=0; i < rows; i++){
      if (w[i] > 0.0){
	wts[i - numzero_weights] = sqrt(w[i]);
	y_wts_f[i - numzero_weights] =  wts[i - numzero_weights]*y[i];
	for (j = 0; j < cols; j++){
	  x_wts_f[j*(rows-totnumzero)+(i-numzero_weights)] = 
	    wts[i - numzero_weights]*x[j*rows+i];
	}
      } else {
	numzero_weights++;
      }
    }
    
    for (j=0;j < cols; j++){
      jpvt[j] = j;
    }
    
    nrows = rows - numzero_weights;
    
    /* now fit the model */
    dqrls_(x_wts_f,&nrows,&cols,y_wts_f,&ny,&tol,beta,resid,qty,&k,jpvt,qraux,work);

    if (k != cols){
      /* the solution is not of full rank */
      /* printf("Improper fit\n"); */
       
      for (j = 0; j < k; j++){
	out_beta[j] = beta[jpvt[j]];
      }
      for(j =k; j < cols; j++){
	out_beta[jpvt[j]] =  R_NaN;
      }

    } else {
    
    /* detangle beta and residual estimates */
      
      for (j = 0; j < cols; j++){
	out_beta[j] = beta[jpvt[j]];
      }
    }
    /* now the model is fitted, lets compute residuals for the 0 weighted observations
       by first computing fitted values. */
    
    numzero_weights = 0;
    for (i=0; i < rows; i++){
      if (w[i] > 0){
	out_resids[i] = resid[i- numzero_weights]/wts[i- numzero_weights];
      } else {
	/* compute the fitted value */
	numzero_weights++;
	fittedvalue = 0.0;
	for (j=0; j <cols; j++){
	  if (out_beta[j] != R_NaN){
	    fittedvalue+=out_beta[j]*x[j*rows+i];
	  }
	}
	out_resids[i] = y[i] -  fittedvalue;
      }
    }

  } else {
    for (i=0; i < rows; i++){
      wts[i] = sqrt(w[i]);
    }
    
    /* dqrls is a linpack routine for least squares solution */
  
    for (i=0; i < rows; i++){
      for (j = 0; j < cols; j++){
	x_wts_f[j*rows+i] = wts[i]*x[j*rows+i];
      }
    }
    
    
    for (i=0; i < rows; i++){
      y_wts_f[i] =  wts[i]*y[i];
    }
    
    for (j=0;j < cols; j++){
      jpvt[j] = j;
    }
    
    /* using function from linpack to fit linear model */


    dqrls_(x_wts_f,&rows,&cols,y_wts_f,&ny,&tol,beta,resid,qty,&k,jpvt,qraux,work);
    
    if (k != cols){
      /* the solution is not of full rank */
      /* printf("Improper fit\n");*/ 
      for (j = 0; j < k; j++){
	out_beta[j] = beta[jpvt[j]];
      }
      for(j =k; j < cols; j++){
	out_beta[j] =  R_NaN; /* out_beta[jpvt[j]] =  R_NaN; */
      }
    } else {
    
      /* detangle beta and residual estimates */
      for (j = 0; j < cols; j++){
	out_beta[j] = beta[jpvt[j]];
      }
    }
    for (i=0; i < rows; i++){
      out_resids[i] = resid[i]/wts[i];
      /* resid[i] = resid[i]/wts[i]; */
    }
  }
  
  Free(wts); 
  Free(x_wts_f); 
  Free(y_wts_f);
  Free(beta);
  Free(resid);
  Free(qraux);
  Free(qty);
  Free(work);
  Free(jpvt );
}


void lvs_rlm_fit(double *x, double *y, int rows, int cols, double *out_beta, double *out_resids, 
	     double *out_weights, int max_iter, int initialized){

  double psi_k = 1.345;
  int i; /* ,j; */
  /* double k = 1.345; */
  /* double k2 = 1.345; */
  double tol = 1e-7;
  double acc = 1e-4;
  double scale =0.0;
  double conv;
  /* int max_iter=20; */
  int iter;
  


  double *wts = out_weights; 
  double *beta = out_beta; 
  double *resids = out_resids; 
  double *old_resids = Calloc(rows,double);
  



  if (!initialized){

    /* intially use equal weights */
    for (i=0; i < rows; i++){
      wts[i] = 1.0;
    }
    
    /* get our intial beta estimates by standard linear regression */
    
    
    lvs_lm_wfit(x, y, wts, rows, cols, tol, beta, resids);
  }

  for (iter = 0; iter < max_iter; iter++){
    
    scale = lvs_med_abs(resids,rows)/0.6745;

    if (fabs(scale) < 1e-10){
      /*printf("Scale too small \n"); */
      break;
    }
    
    for (i =0; i < rows; i++){
      old_resids[i] = resids[i];
    }

    for (i=0; i < rows; i++){
      wts[i] = lvs_psi_huber(resids[i]/scale,psi_k,0);
    }
   
    lvs_lm_wfit(x, y, wts, rows, cols, tol, beta, resids);


    /*check convergence  based on residuals */
    
    conv = lvs_irls_delta(old_resids,resids, rows);

    if (conv < acc){
      /*    printf("Converged \n");*/
      break; 
      
    }
  }


  Free(old_resids);
}


void lvs_rlm_fit_R(double *x, double *y, int *rows, int *cols, double *out_beta, 
		   double *out_resids, double *out_weights)
{
  lvs_rlm_fit(x, y, *rows, *cols, out_beta, out_resids,out_weights,20,0);
}

