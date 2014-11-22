/*
 * initest.h
 *
 *  Created on: Jul 13, 2014
 *      Author: Rafat Hussain
 */

#ifndef INITEST_H_
#define INITEST_H_

#include "boxjenkins.h"

#ifdef __cplusplus
extern "C" {
#endif

void innalg(double *vec, int l, int N,double *v, double *theta);

void ma_inn(double *inp, int N,int q, double* theta, int lag);

void burgalg(double *x, int N, int p, double *phi,double *var);

void ywalg(double *x, int N, int p, double *phi);

void ywalg2(double *x, int N, int p, double *phi, double *var);

void hralg(double *x, int N, int p,int q, double *phi,double *theta, double *var);

void hrstep2(double *inp, int N, int m, int p, int q, int pq, double *a, double *sos, double *var);

void pacf_burg(double* vec, int N, double* par, int M);

void pacf_yw(double* vec, int N, double* par, int M);

void pacf_mle(double* vec, int N, double* par, int M);

#ifdef __cplusplus
}
#endif


#endif /* INITEST_H_ */
