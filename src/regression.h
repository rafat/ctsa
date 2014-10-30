/*
 * regression.h
 *
 *  Created on: Jun 5, 2013
 *      Author: USER
 */

#ifndef REGRESSION_H_
#define REGRESSION_H_

#include "stats.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct bparam_t {
  double value;
  double lower;
  double upper;
  double stdErr;
} bparam;

typedef struct reg_set* reg_object;

reg_object reg_init(int N, int p);

struct reg_set{
	int N;
	int p;
	double alpha;
	double sigma;
	double sigma_lower;
	double sigma_upper;
	double R2[2];
	int df;
	double TSS;
	double ESS;
	double RSS;
	int df_ESS;
	int df_RSS;
	double FStat;
	double PVal;
	bparam beta[1];
};

void linreg_clrm(double *x,double *y, int N, double* b,
		double *var,double *res,double alpha,double *anv,
		double* ci_lower, double* ci_upper);
		
void zerohyp_clrm(int N,double *b, double *val, double *tval, double *pval);

//void linreg_multi2(int p, double *x,double *y, int N, double* b); // p number of variables.
// p = 2 for one dependent variable	and one independent variable
// p = 3 for one dependent variable	and two independent variables etc.

void linreg_multi(int p, double *x,double *y, int N, double* b,double *sigma2,
			double *xxti,double *R2,double *res,double alpha,double *anv,
		double* ci_lower, double* ci_upper);
		
void zerohyp_multi(int N,double *b,int p, double *varcovar, double *tval, double *pval);

void regress(reg_object obj,double *x,double *y,double *res,double *varcovar,double alpha);

void regress_poly(reg_object obj,double *x,double *y,double *res,double *varcovar,double alpha);

void summary(reg_object obj);

void anova(reg_object obj);

void confint(reg_object obj);

void zerohyp_val(reg_object obj, double *tval, double *pval);

double fitted(reg_object obj,double *inp,double *varcovar,double *var);

void free_reg(reg_object object);		


#ifdef __cplusplus
}
#endif

#endif /* REGRESSION_H_ */
