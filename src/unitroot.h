#ifndef UNITROOT_H_
#define UNITROOT_H_

#include "errors.h"

#ifdef __cplusplus
extern "C" {
#endif

void ur_df(double *y, int N,const char* alternative, int *klag, double *statistic,double *pval);

void ur_kpss(double *y, int N,const char* type,int lshort, int *klag, double *statistic,double *pval);

void ur_pp(double *y, int N,const char* alternative,const char* type,int lshort, int *klag, double *statistic,double *pval);

void ur_pp2(double *x, int N,const char* type,const char* model,int lshort, int *klag,double *cval,double *cprobs, double *auxstat,int *laux,double *teststat);

void ur_df2(double *y, int N,const char* type, int *lags,const char *selectlags,double *cval,int *cvrows, int *cvcols, double *cprobs, double *teststat,int *ltstat);

#ifdef __cplusplus
}
#endif


#endif /* UNITROOT_H_ */