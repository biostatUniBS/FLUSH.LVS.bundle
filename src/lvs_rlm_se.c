/*Time-stamp: <04-06-2009 19:55:32 ste DrSlump> */

/*********************************************************************
 **
 ** file: lvs_rlm_se.c
 **
 ** based on rlm_se.c from Ben Bolstad's affyPLM library (version 1.4.0)
 **
 ** Aim: implement computation of standard errors for robust linear models.
 **
 ** Copyright (C) 2003-2004 Ben Bolstad
 **
 ** created on: May 10, 2007
 **
 ** Last modified: June 4, 2009
 **
 **  
 ********************************************************************/

#include "basic_fns.h"
#include "lvs_matrix_functions.h"
#include "externals.h"
#include "lvs_rlm_se.h"

#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>





/*********************************************************************
 **
 ** void lvs_rlm_compute_se(double *X,double *Y, int n, int p, double *beta, double *resids,double *weights,double *se_estimates, int method)
 **
 ** given a robust linear model fit to data, compute the standard errors of the parameter estimates
 **
 ** double *X
 ** double *Y
 ** int n
 ** int p
 ** double *beta
 ** double *resids
 ** double *weights
 ** double *se_estimates
 ** int method
 **
 **
 ** Formulas refer to Huber 2nd Ed.
 **
 **
 ** Note that we compute Kappa using Huber's eq 7.84 pag 171
 **
 ** Kappa = 1 + p/n* var(psi')/(E psi')^2 
 **
 **  
 ** 
 **
 **
 *********************************************************************/

void lvs_rlm_compute_se(double *X,double *Y, int n, int p, double *beta, double *resids,
			double *weights,
			double *residSE,double *out_X,
			int start_p, int n_b)
{
  int i,j,k;

  /* double resvar = 0.0; */
  double scale = 0.0;
  double Kappa = 0.0;
  double sumderivpsi=0.0; /* sum of psi'(r_i) */
  double sumpsi2 = 0.0;  /* sum of psi(r_i)^2 */
  double vs=0.0,m,varderivpsi=0.0;
  double RMSEw = 0.0;

  double k_huber = 1.345;

  double *XTX = Calloc(p*p,double);
  
  scale = lvs_med_abs(resids,n)/0.6745;
  
/*   for (i=0; i < n; i++) */
/*     { */
/*       RMSEw+= weights[i]*resids[i]*resids[i]; */
/*     } */
  
/*   RMSEw = sqrt(RMSEw/(double)(n-p)); */


  for (i =0; i < n; i++){
    //sum psi(r_i)^2
    sumpsi2+= lvs_psi_huber(resids[i]/scale,k_huber,2)*lvs_psi_huber(resids[i]/scale,k_huber,2);
    // sum psi'(r_i)
    sumderivpsi+= lvs_psi_huber(resids[i]/scale,k_huber,1);
  }
  
  m = (sumderivpsi/(double) n); // E[psi(r)']
  
  
  for (j =0; j < p;j++)
    {
      for (k=0; k < p; k++){
	XTX[k*p + j] = 0.0;
	for (i=0; i < n; i++){
	  XTX[k*p + j]+=  X[j*n +i]*X[k*n + i];
	}
      }

    }
  
  for (i = 0; i < n; i++){
    varderivpsi+=(lvs_psi_huber(resids[i]/scale,k_huber,1) - m)
      *(lvs_psi_huber(resids[i]/scale,k_huber,1) - m);
  }

  varderivpsi/=(double)(n);
  

  Kappa = 1.0 + (double)p *varderivpsi/((double)n*m*m);

  vs = scale*scale*sumpsi2/(double)(n-p);
  
  residSE[0] =  Kappa*sqrt(vs)/m;

  Kappa = residSE[0]*residSE[0];
  
  /* residSE[0] = RMSEw; */
  /* Kappa = RMSEw*RMSEw; */
  
  estimate_X(beta, start_p, n_b, p, XTX, Kappa,out_X);
  
  Free(XTX);
}


void estimate_X(double *beta, int start_p, int n_b, int p, double *XTX, double Kappa,
		double *out_X)
{
  int i = 0;
  int j = 0;
  int end_p = 0.0;
  
  out_X[0] = 0.0;

  /* to avoid thinking in terms of [0] first pointer in R code*/
  start_p = start_p - 1; /* position first beta sample*/
  end_p = start_p + n_b ; /* position last beta sample. n_b = number betas for sample*/

  
  for(i = start_p; i < end_p ; i++)
    {
      out_X[0] += (beta[i]*beta[i]* XTX[i*p + i]) / Kappa;
    }
}



/* R wrapper functions  */




void lvs_rlm_compute_se_R(double *X,double *Y, int *n, int *p, double *beta, double *resids,
			  double *weights, double *residSE,
			  double *out_X,  int *start_p,int *n_b)
{
  lvs_rlm_compute_se(X, Y, *n, *p, beta, resids, weights, residSE,out_X, *start_p, *n_b);
  
}



void estimate_X_R(double *beta, int *start_p, int *n_b, int *p, double *XTX, double *Kappa,
		  double *out_X)
{
  estimate_X(beta, *start_p, *n_b, *p, XTX, *Kappa,out_X);
}
