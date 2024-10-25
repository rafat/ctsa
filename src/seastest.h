// SPDX-License-Identifier: BSD-3-Clause
#ifndef SEASTEST_H_
#define SEASTEST_H_

#include "unitroot.h"

#ifdef __cplusplus
extern "C" {
#endif

void decompose(double *x,int N, int f,double *filter, const char *type, double *trend, int *ltrend, double *seas, int *lseas, double *random, int *lrandom);

double* genLags(double *y, int N, int maxLags, int *rows, int *cols);

reg_object fitOCSB(double *x, int N, int f, int lag, int mlags);

void OCSBtest(double *x, int N, int f, int mlags, const char *method,double *statistics,double *critical);

void stl(double *x,int N,int f, const char *s_window_type,int *s_window, int *s_degree, int *t_window, int *t_degree,int *l_window,int *l_degree,
    int *s_jump, int *t_jump, int *l_jump, int *robust,int *inner, int *outer,double *seasonal,double *trend, double *remainder);

void modstl(double *x, int N, int f, int *s_window,double *lambda, double *seasonal, double *trend,double *remainder);

void mstl(double *x, int N, int *f, int *Nseas, int *s_window,double *lambda,int *iterate, double **seasonal, double *trend,double *remainder);

void SHtest(double *x, int N, int *f, int Nseas, double *season);

double* seasdummy(double *x, int N,int f,int *rows, int *cols);

#ifdef __cplusplus
}
#endif


#endif /* SEASTEST_H_ */