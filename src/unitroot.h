#ifndef UNITROOT_H_
#define UNITROOT_H_

#include "errors.h"

#ifdef __cplusplus
extern "C" {
#endif

void ur_df(double *y, int N,const char* alternative, int *klag, double *statistic,double *pval);

void ur_kpss(double *y, int N,const char* type,int lshort,double *statistic,double *pval);

#ifdef __cplusplus
}
#endif


#endif /* UNITROOT_H_ */