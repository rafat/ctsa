#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../header/ctsa.h"

int main(void) {
	int i, N, L,method;
	double *inp;
	int p;
	double *phi;
	double *xpred, *amse;

	ar_object obj;
	p = 0;


	L = 5;

	phi = (double*)malloc(sizeof(double)* p);

	xpred = (double*)malloc(sizeof(double)* L);
	amse = (double*)malloc(sizeof(double)* L);

	FILE *ifp;
	double temp[2000];

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
	//wmean = mean(temp, N);

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
		//printf("%g \n",inp[i]);
	}

	method = 0; // method 0 - Yule Walker, Method 1 - Burg, Method 2, MLE (Box-Jenkins)
	obj = ar_init(method, N);
	ar_exec(obj, inp);
	ar_summary(obj);
	// Predict the next 5 values using the obtained ARIMA model
	ar_predict(obj, inp, L, xpred, amse);
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

	ar_free(obj);
	//ar_estimate(inp, N, 2);
	free(inp);
	free(phi);
	free(xpred);
	free(amse);
	return 0;

}
