#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../header/ctsa.h"
/*

ARMA example for a stationary time series with no seasonal components and no exogenous variable.

*/
int main(void) {
	int i, N, d, D, L;
	double *inp;
	int p, q, P, Q, s, r,imean;
	double *phi, *theta;
	double *xpred, *amse;
	sarimax_object obj;
    /*
    Make sure all the parameter values are correct and consistent with other values. eg., if xreg is NULL r should be 0
    or if P = D = Q = 0 then make sure that s is also 0. 
     Recheck the values if the program fails to execute.
    */
	p = 0;
	d = 1;
	q = 0;
    s = 0;
	P = 0;
	D = 0;
	Q = 0;
	r = 0;

	imean = 1;

	L = 5;

	phi = (double*)malloc(sizeof(double)* p);
	theta = (double*)malloc(sizeof(double)* q);

	xpred = (double*)malloc(sizeof(double)* L);
	amse = (double*)malloc(sizeof(double)* L);

	FILE *ifp;
	double temp[1200];

	ifp = fopen("../data/seriesA.txt", "r");
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

	inp = (double*)malloc(sizeof(double)* N);

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
	}


	obj = sarimax_init(p, d, q, P, D, Q, s, r ,imean, N);

    /* setMethod()
    Method 0 ("CSS-MLE") is default. The method also accepts values 1 ("MLE") and 2 ("CSS")
    */

	sarimax_setMethod(obj, 0); 

    /*sarimax_exec(object, input time series, exogenous time series)
        set exogenous to NULL if deadling only with a univariate time series.
    */
	sarimax_exec(obj, inp,NULL);
	sarimax_summary(obj);
	/* sarimax_predict(sarimax_object obj, double *inp, double *xreg, int L,double *newxreg, double *xpred, double *amse)
        inp - Input Time Series
        xreg - Exogenous Time Series
        L - L point prediction
        newxreg - Exogenous Time Series of dimension r * L where r is the number of exogenus time series and L is the length of each
        xpred - L future values
        amse - MSE for L future values
    */
	sarimax_predict(obj, inp, NULL, L, NULL, xpred, amse);
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
	free(phi);
	free(theta);
	free(xpred);
	free(amse);
	return 0;

}
