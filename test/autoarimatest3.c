#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../header/ctsa.h"
/*

Auto ARIMA example 1

*/

int main(void) {
    int i, N, d, D, L;
	double *inp;
	int p, q, P, Q, s, r;
	double *xpred, *amse,*xreg,*newxreg;
	auto_arima_object obj;
    /*
    Make sure all the parameter values are correct and consistent with other values. eg., if xreg is NULL r should be 0
    or if P = D = Q = 0 then make sure that s is also 0. 
     Recheck the values if the program fails to execute.
    */
	p = 5;
	d = 2;
	q = 5;
	P = 2;
	D = 1;
	Q = 2;
	
	int order[3] = {p,d,q};
	int seasonal[3] = {P,D,Q};

    s = 0;
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


	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
        xreg[i] = temp1[i];
		xreg[N+i] = temp2[i];
	}

    for(i = 0; i < L;++i) {
        newxreg[i] = temp1[N + i];
        newxreg[i+L] = temp2[N + i];
    }

	obj = auto_arima_init(order,seasonal,s,r,N);

	auto_arima_setApproximation(obj,0);
	auto_arima_setStepwise(obj,0);
	auto_arima_setVerbose(obj,1);

	auto_arima_exec(obj,inp,xreg);

	auto_arima_summary(obj);

	auto_arima_predict(obj, inp, xreg, L,newxreg, xpred, amse);

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

	auto_arima_free(obj);

	
	free(inp);
	free(xpred);
	free(amse);
    free(xreg);
    free(newxreg);

    return 0;
}