#ifndef BOXCOX_H_
#define BOXCOX_H_

#include "optimc.h"
#include "stats.h"

#ifdef __cplusplus
extern "C" {
#endif

int checkConstant(double *x, int N);

void boxcox_eval(double *x, int N, double lambda,double *bxcx);

double inv_boxcox_eval(double *x,int N, double lambda,double *bxcx);

double boxcox_loglik(double lambda, void *params);

void boxcox(double *x, int N,double *lambda,double *y);

#ifdef __cplusplus
}
#endif



#endif /* BOXCOX_H_ */