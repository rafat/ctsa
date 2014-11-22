/*
 * boxjenkins.h
 *
 *  Created on: Jun 29, 2014
 *      Author: Rafat Hussain
 */

#ifndef BOXJENKINS_H_
#define BOXJENKINS_H_

#include "talg.h"
#include "optimc.h"


#ifdef __cplusplus
extern "C" {
#endif

void USPE(double *inp,int N,int p, int q, double *phi,double *theta,double *thetac,double *var);

void USPE_seasonal(double *inp,int N,int s,int ps, int qs, double *phis,double *thetas);

void avaluem(double *w,int N,int p, int q, double *phi,double *theta,int tval,double *a) ;

int nlalsm(double *vec, int N, int p, int delta, int q, double *phi, double *theta,
	int M, double *thetac, double *var, double eps, double *varcovar, double *residuals);

int nlalsms(double *vec,int N,int p,int delta, int q, double *phi,double *theta,
		int s,int lps,int deltas,int lqs,double *phis,double *thetas,
		int M,double *thetac,double *var,double eps,double *varcovar,double *residuals);


#ifdef __cplusplus
}
#endif


#endif /* BOXJENKINS_H_ */
