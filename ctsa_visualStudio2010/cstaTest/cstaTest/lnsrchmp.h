#ifndef LNSRCHMP_H_
#define LNSRCHMP_H_

#include "matrix.h"

#define EPSILON 2.7182818284590452353602874713526624977572

#ifdef __cplusplus
extern "C" {
#endif

typedef struct custom_function_set custom_function;

struct custom_function_set{
	double(*funcpt) (double *x,int N,void *params);// Function in N variables
	void *params;
};

typedef struct custom_gradient_set custom_gradient;

struct custom_gradient_set{
	void(*funcgrad) (double *x, int N,double *g, void *params);
	void *params;
};

typedef struct custom_funcuni_set custom_funcuni;

struct custom_funcuni_set{
	double(*funcuni) (double x, void *params);// Function in one variable
	void *params;
};

typedef struct custom_funcmult_set custom_funcmult;

struct custom_funcmult_set{
	void(*funcmult) (double *x,int M, int N, double *f, void *params);// M Functions in N variables
	void *params;
};

typedef struct custom_jacobian_set custom_jacobian;

struct custom_jacobian_set{
	void(*jacobian) (double *x, int M, int N, double *jac, void *params);
	void *params;
};

#define FUNCPT_EVAL(F,x,N) (*((F)->funcpt))(x,N,(F)->params)

#define FUNCGRAD_EVAL(F,x,N,g) (*((F)->funcgrad))(x,N,(g),(F)->params)

#define FUNCUNI_EVAL(F,x) (*((F)->funcuni))(x,(F)->params)

#define FUNCMULT_EVAL(F,x,M,N,f) (*((F)->funcmult))(x,M,N,(f),(F)->params)

#define JACOBIAN_EVAL(F,x,M,N,jac) (*((F)->jacobian))(x,M,N,(jac),(F)->params)

int stopcheck_mt(double fx, int N, double *xc, double *xf, double *jac, double *dx, double fsval, double gtol, double stol, int retval);

int stopcheck2_mt(double fx, int N, double fo, double *jac, double *dx, double eps,double stoptol, double functol, int retval);

int stopcheck3_mt(double *xi,double *xf,double fx, int N, double fo, double *jac, double *dx, double eps,
		double stoptol, double functol, int retval);

int grad_fd(custom_function *funcpt,custom_gradient *funcgrad, double *x, int N, double *dx, double eps2, double *f);

int grad_cd(custom_function *funcpt, custom_gradient *funcgrad, double *x, int N, double *dx,
		double eps3, double *f);

int grad_calc2(custom_function *funcpt, double *x, int N, double *dx, double eps3, double *f);

int grad_calc(custom_function *funcpt, double *x, int N, double *dx, double eps2, double *f);

int cstep(double *stx, double *fx, double *dx, double *sty, double *fy, double *dy, double *stp, double *fp, double *dp, int *brackt,
	double  stpmin, double stpmax);

int cvsrch(custom_function *funcpt, custom_gradient *funcgrad, double *x, double *f, double *g, double *stp, double *s, int N, double *dx, double maxstep,
	int MAXITER,double eps2,double ftol, double gtol, double xtol);

int lnsrchmt(custom_function *funcpt, custom_gradient *funcgrad, double *xi, double *f, double *jac, double *alpha, double *p, int N, double *dx, double maxstep, int MAXITER,
		double eps2,double ftol, double gtol, double xtol, double *x);

#ifdef __cplusplus
}
#endif

#endif /* LNSRCHMP_H_ */
