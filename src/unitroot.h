#ifndef UNITROOT_H_
#define UNITROOT_H_

#include "errors.h"

#ifdef __cplusplus
extern "C" {
#endif

void ur_df(double *y, int N,const char* type, int lags, const char *selectlags,double *statistic);

#ifdef __cplusplus
}
#endif


#endif /* UNITROOT_H_ */