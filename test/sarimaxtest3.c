// SPDX-License-Identifier: BSD-3-Clause
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../header/ctsa.h"
/*

ARIMA example with exogenous variables


*/
int main(void) {
	int i, N, d, D, L;
	double *inp;
	int p, q, P, Q, s, r;
	double *xpred, *amse,*xreg,*newxreg;
	sarimax_object obj;
	int imean = 1;
    /*
    Make sure all the parameter values are correct and consistent with other values. eg., if xreg is NULL r should be 0
    or if P = D = Q = 0 then make sure that s is also 0. 
     Recheck the values if the program fails to execute.
    */
	p = 2;
	d = 0;
	q = 2;
	s = 0;
	P = 0;
	D = 0;
	Q = 0;
	r = 2;


	L = 5;

	xpred = (double*)malloc(sizeof(double)* L);
	amse = (double*)malloc(sizeof(double)* L);

	FILE *ifp;
	double temp[1200];
    double temp1[1200];
    double temp2[1200];

	ifp = fopen("../data/e1m.dat", "r");
	i = 0;
	if (!ifp) {
		printf("Cannot Open File");
		exit(100);
	}
	while (!feof(ifp)) {
		fscanf(ifp, "%lf %lf %lf \n", &temp[i],&temp1[i],&temp2[i]);
		i++;
	}
	N = i - L;

	inp = (double*)malloc(sizeof(double)* N);
    xreg = (double*)malloc(sizeof(double)* N * 2);
    newxreg = (double*)malloc(sizeof(double)* L * 2);

    /*
    
    */

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
        xreg[i] = temp1[i];
		xreg[N+i] = temp2[i];
	}

    for(i = 0; i < L;++i) {
        newxreg[i] = temp1[N + i];
        newxreg[i+L] = temp2[N + i];
    }

	obj = sarimax_init(p, d, q, P, D, Q, s, r ,imean, N);

    /* setMethod()
    Method 0 ("CSS-MLE") is default. The method also accepts values 1 ("MLE") and 2 ("CSS")
    */

	sarimax_setMethod(obj, 2); 

    /*sarimax_exec(object, input time series, exogenous time series)
        set exogenous to NULL if deadling only with a univariate time series.
    */
	sarimax_exec(obj, inp,xreg);
	sarimax_summary(obj);
	/* sarimax_predict(sarimax_object obj, double *inp, double *xreg, int L,double *newxreg, double *xpred, double *amse)
        inp - Input Time Series
        xreg - Exogenous Time Series
        L - L point prediction
        newxreg - Exogenous Time Series of dimension r * L where r is the number of exogenus time series and L is the length of each
        xpred - L future values
        amse - MSE for L future values
    */
	sarimax_predict(obj, inp, xreg, L, newxreg, xpred, amse);
	printf("\n");
	printf("Predicted Values : ");
	for (i = 0; i < L; ++i) {
		printf("%g ", xpred[i]);
	}
	printf("\n");
	printf("Standard Errors  : ");
	for (i = 0; i < L; ++i) {
		printf("%g ", sqrt(amse[i]));
	}
	printf("\n");

	sarimax_free(obj);
	free(inp);
	free(xpred);
	free(amse);
    free(xreg);
    free(newxreg);
	return 0;

}
