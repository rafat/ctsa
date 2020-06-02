#ifndef AUTOUTILS_H_
#define AUTOUTILS_H_

#include "seastest.h"

#ifdef __cplusplus
extern "C" {
#endif

int ndiffs(double *x, int N,double *alpha, const char *test, int *max_d);

int nsdiffs(double *x, int N,int f,double *alpha, const char *test, int *max_D);


#ifdef __cplusplus
}
#endif


#endif /* AUTOUTILS_H_ */