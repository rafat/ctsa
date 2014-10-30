/*
 * nls.h
 *
 *  Created on: May 21, 2014
 *      Author: HOME
 */

#ifndef NLS_H_
#define NLS_H_

#include "lls.h"

#ifdef __cplusplus
extern "C" {
#endif

double enorm(double *x, int N);

void qrfac(double *A, int M, int N, int lda, int pivot, int *ipvt, int lipvt,double *rdiag, double *acnorm,double eps);

void qrsolv(double *r,int ldr,int N,int *ipvt,double *diag,double *qtb,double *x,double *sdiag);

void fdjac2(custom_funcmult *funcmult, double *x, int M, int N, double *fvec, double *fjac, int ldfjac,
		double epsfcn,double eps);

void lmpar(double *r,int ldr,int N,int *ipvt,double *diag,double *qtb,double delta,double *par,double *x,double *sdiag);

int lmder(custom_funcmult *funcmult,custom_jacobian *jacobian,double *xi,int M, int N,
		double *fvec,double *fjac,int ldfjac,int maxfev,double *diag,int mode,double factor,int nprint,
		double eps,double ftol,double gtol,double xtol,int *nfev,int *njev,int *ipvt, double *qtf);

int lmdif(custom_funcmult *funcmult, double *x, int M, int N, double *fvec, double *fjac, int ldfjac,
		int maxfev,double *diag,int mode,double factor,int nprint,double eps,double epsfcn,double ftol,double gtol,
		double xtol,int *nfev,int *njev,int *ipvt, double *qtf);


#ifdef __cplusplus
}
#endif


#endif /* NLS_H_ */
