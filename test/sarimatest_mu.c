// SPDX-License-Identifier: BSD-3-Clause
/**
 * Adapting to minunit test library
 * Started on Fri Jan 17 11:57:23 2025 
 * From: `sarimatest.c` 
 **/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../header/ctsa.h"

#include "minunit.h" /* Test support macros */

struct expected_results {
	double aic;
	double log_likelihood;
	double sigma_square;
};

struct singleton_like {
	double tolerance;
	double *inp;
	double *xpred;
	double *amse;
	sarima_object obj;
	struct expected_results expected;	
};

struct singleton_like globals = {0};

void test_setup(void) {
	int i;
	int N;	
	int p;
	int q;
	int d;	
	int L;	
	int P;	
	int D;	
	int Q;	
	int s;	
	
	FILE *ifp;
	double temp[1200];
	
	p = 0;
	d = 1;
	q = 1;
	P = 0;
	D = 1;
	Q = 1;
	s = 12;
	L = 5;
	
	globals.xpred = (double*)(double*)calloc(sizeof(double), L);
	globals.amse = (double*)(double*)calloc(sizeof(double), L);

	ifp = fopen("../data/seriesG.txt", "r");
	i = 0;
	if (!ifp) {
		printf("Cannot Open File");
		exit(100);
	}
	while (!feof(ifp)) {
		fscanf(ifp, "%lf \n", &temp[i]);
		i++;
	}
	N = i;
	
	fclose(ifp);

	globals.inp = (double*)calloc(sizeof(double), N);

	for (i = 0; i < N; ++i) {
		globals.inp[i] = log(temp[i]);
		//printf("%g \n",inp[i]);
	}

	globals.obj = sarima_init(p, d, q, s, P, D, Q, N);
	
	sarima_setMethod(globals.obj, 0); 
	/** 
	 * Method 0 ("MLE") is default so this step is unnecessary. The 
	 * 	method also accepts values 1 ("CSS") and 2 ("Box-Jenkins") 
	 **/
	
	//sarima_setOptMethod(globals.obj, 7);
	/** 
	 * Method 7 ("BFGS with More Thuente Line Search") is default so 
	 * this step is unnecessary. The method also accepts values 
	 * 0,1,2,3,4,5,6. Check the documentation for details.
	 **/
	 
	 /** Execute the fitting of the SARIMA model */
	sarima_exec(globals.obj, globals.inp);
	 
	 /** Set tolerance */	
	globals.tolerance = 0.01;
	 
	 /** Expected results assignment */
	globals.expected.aic = -483.145;
	globals.expected.log_likelihood = 244.572;
	globals.expected.sigma_square = 0.00135071;
}

void test_teardown(void) {
	sarima_free(globals.obj);
	
	free(globals.inp);
	free(globals.xpred);
	free(globals.amse);
}

MU_TEST(test_assert_var) {	
	sarima_summary(globals.obj);
	
	mu_assert_double_tol(globals.expected.sigma_square, 
						 globals.obj->var, 
						 globals.tolerance);
}

MU_TEST(test_assert_loglik) {
	mu_assert_double_tol(globals.expected.log_likelihood, 
						 globals.obj->loglik, 
						 globals.tolerance);	
}

MU_TEST(test_assert_aic) {
	mu_assert_double_eq(globals.expected.aic, globals.obj->aic);
}

MU_TEST_SUITE(test_suite) {
	MU_SUITE_CONFIGURE(&test_setup, &test_teardown);
	
	MU_RUN_TEST(test_assert_var);
	MU_RUN_TEST(test_assert_loglik);
	MU_RUN_TEST(test_assert_aic);
}

int main(void) {
	MU_RUN_SUITE(test_suite);
	MU_REPORT();
	return MU_EXIT_CODE;
}
