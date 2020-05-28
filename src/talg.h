/*
 * talg.h
 *
 *  Created on: Jul 31, 2013
 *      Author: Rafat Hussain
 */

#ifndef TALG_H_
#define TALG_H_

#include "filter.h"
#include "regression.h"
#include "conv.h"
#include "polyroot.h"

#ifdef __cplusplus
extern "C" {
#endif

void detrend_ma(double *sig,int N, int window,double *oup);

int poly(double *A, double *B, double *C, int lA, int lB);

int upsample(double *x,int lenx, int M, double *y);

int downsample(double *x,int lenx, int M, double *y);

void deld(int d, double* oup);

void delds(int D, int s, double *C);

int diff(double *sig, int N, int del, double *oup);

int diffs(double *sig, int N, int D,int s, double *oup);

void deseason_ma(double *sig,int N,int s,double *oup);

void psiweight(double *phi,double *theta,double *psi,int p,int q,int j);

void piweight(double *phi,double *theta,double *piw,int p,int q,int j);

void arma_autocovar(double *phi,double *theta,int p,int q,double var,double* acov, int lag);

int twacf(double *P, int MP, double *Q, int MQ, double *ACF, int MA, double *CVLI, int MXPQ1, double *ALPHA, int MXPQ);

void artrans(int p, double *old, double *new1);

void arinvtrans(int p, double *old, double *new1);

void transall(int p, int q, int P, int Q, double *old, double *new1);

void invtransall(int p, int q, int P, int Q, double *old, double *new1);

int invertroot(int q, double *ma);

double interpolate_linear(double *x,double *y, int N, double z);

void approx(double *x,double *y, int N,double *xout, double *yout,int Nout);

void linspace(double *x, int N,double xlo,double xhi);

void arrayminmax(double *x, int N, double *amin,double *amax);

void cumsum(double *x, int N, double *csum);

void ppsum(double *x, int N, int l, double *sum);

void stl_(double *y, int *n, int *np, int *ns, int *nt, int *nl, int *isdeg, int *itdeg, int *
	ildeg, int *nsjump, int *ntjump, int *nljump, int *ni, int *no, double *rw, double *season, double *trend, 
	double *work);

int stless_(double *y, int *n, int *len, int *ideg, int *njump, int *userw, double *rw, double *ys,
	 double *res);

int stlest_(double *y, int *n, int *len, int *ideg, double *xs, double *ys, int *nleft, int *
	nright, double *w, int *userw, double *rw, int *ok);

int stlfts_(double *x, int *n, int *np, double *trend, double *work);

int stlma_(double *x, int *n, int *len, double *ave);

int stlstp_(double *y, int *n, int *np, int *ns, int *nt, int *nl, int *isdeg, int *itdeg, int 
	*ildeg, int *nsjump, int *ntjump, int *nljump, int *ni, int *userw, double *rw, double *season,
	double *trend, double *work);

int stlrwt_(double *y, int *n, double *fit, double *rw);

int stlss_(double *y, int *n, int *np, int *ns, int *isdeg, int *nsjump, int *userw, double *rw, 
	double *season, double *work1, double *work2, double *work3, double *work4);

void stlez_(double *y, int *n, int *np, int *ns, int *isdeg, int *itdeg,int *robust, int *no, 
	double *rw, double *season, double *trend, double *work);

int psort_(double *a, int n, int *ind, int ni);

#ifdef __cplusplus
}
#endif



#endif /* TALG_H_ */
