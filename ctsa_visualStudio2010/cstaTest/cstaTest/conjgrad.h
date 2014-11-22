#ifndef CONJGRAD_H_
#define CONJGRAD_H_

#include "newtonmin.h"


#ifdef __cplusplus
extern "C" {
#endif

int ichol(double *A, int N);

int stopcheck2(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol) ;

int cgpr_mt(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx,double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double ftol,double xtol,double *xf); //Polak Ribiere + (More Thuentes Line Search)

int conjgrad_min_lin(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double ftol,double xtol,double *xf);


#ifdef __cplusplus
}
#endif

#endif /* CONJGRAD_H_ */
