// cstaTest.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ctsa.h"





int _tmain(int argc, _TCHAR* argv[])
{
	
	int i, N, d, L;
	double *inp;
	int p, q;
	double *phi, *theta;
	double *xpred, *amse;
	arima_object obj;
	p = 1;
	d = 0;
	q = 1;


	L = 5;

	phi = (double*)malloc(sizeof(double)* p);
	theta = (double*)malloc(sizeof(double)* q);

	xpred = (double*)malloc(sizeof(double)* L);
	amse = (double*)malloc(sizeof(double)* L);

	FILE *ifp;
	double temp[1200];

	ifp = fopen("seriesA.txt", "r");
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
	//wmean = mean(temp, N);

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
		//printf("%g \n",inp[i]);
	}


	obj = arima_init(p, d, q, N);
	arima_setMethod(obj, 0); // Method 0 ("MLE") is default so this step is unnecessary. The method also accepts values 1 ("CSS") and 2 ("Box-Jenkins")
	arima_setOptMethod(obj, 7);// Method 7 ("BFGS with More Thuente Line search") is default so this step is unnecessary. The method also accepts values 0,1,2,3,4,5,6. Check the documentation for details.
	arima_exec(obj, inp);
	arima_summary(obj);
	// Predict the next 5 values using the obtained ARIMA model
	arima_predict(obj, inp, L, xpred, amse);
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
	arima_free(obj);
	free(inp);
	free(phi);
	free(theta);
	free(xpred);
	free(amse);
	
	
	
	
	
	return 0;
}

