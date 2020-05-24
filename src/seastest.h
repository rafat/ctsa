#ifndef SEASTEST_H_
#define SEASTEST_H_

#include "unitroot.h"

#ifdef __cplusplus
extern "C" {
#endif

void decompose(double *x,int N, int f,double *filter, const char *type, double *trend, int *ltrend, double *seas, int *lseas, double *random, int *lrandom);

double* genLags(double *y, int N, int maxLags, int *rows, int *cols);

reg_object fitOCSB(double *x, int N, int f, int lag, int mlags);

void OCSBtest(double *x, int N, int f, int mlags, const char *method);


#ifdef __cplusplus
}
#endif


#endif /* SEASTEST_H_ */