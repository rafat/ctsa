#ifndef ERRORS_H_
#define ERRORS_H_

#include "initest.h"

#ifdef __cplusplus
extern "C" {
#endif

double me(double *predicted, double *actual, int N);

double mse(double *predicted, double *actual, int N);

double mae(double *predicted, double *actual, int N);

double mape(double *predicted, double *actual, int N);

double mpe(double *predicted, double *actual, int N);

double mase(double *predicted, double *actual, int N, double *tseries, int length);

#ifdef __cplusplus
}
#endif


#endif /* ERRORS_H_ */
