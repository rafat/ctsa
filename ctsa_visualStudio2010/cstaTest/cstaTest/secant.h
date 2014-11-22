#ifndef SECANT_H_
#define SECANT_H_

#include "conjgrad.h"


#ifdef __cplusplus
extern "C" {
#endif

void bfgs_naive(double *H,int N,double eps,double *xi,double *xf,double *jac,double *jacf);

int bfgs_min_naive(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double fsval,double maxstep, int MAXITER,
		double eps, double *xf);

void bfgs_factored(double *H,int N,double eps,double *xi,double *xf,double *jac,double *jacf);

int bfgs_min(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double fsval, double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double stol,double *xf);

int bfgs_min2(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, int m, double *dx, double fsval, double maxstep, int MAXITER, int *niter,
	double eps, double gtol, double ftol, double xtol, double *xf);

void inithess_l(double *H, int N, int k, double *tsk, double *tyk, double *dx);

void bfgs_rec(double *H, int N, int iter, int m, double *jac, double *sk, double *yk, double *r);

int bfgs_l_min(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, int m, double *dx, double fsval,double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double ftol,double xtol,double *xf);

#ifdef __cplusplus
}
#endif

#endif /* SECANT_H_ */
