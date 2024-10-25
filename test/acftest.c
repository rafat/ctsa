// SPDX-License-Identifier: BSD-3-Clause
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../header/ctsa.h"

int main(void) {
	int i, N, M;
	double *inp, *acf;
	int method;

	/*
	acvf and acvf_opt functions will calculate autocovariance from lag 0 to lag M-1
	*/

	FILE *ifp;
	double temp[1200];

	ifp = fopen("../data/seriesC.txt", "r");
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

	M = 10;
	acf = (double*)malloc(sizeof(double)* M);

	// Default Method

	printf("\n Default Method : acvf \n");

	acvf(inp, N, acf, M);

	for (i = 0; i < M; ++i) {
		printf("%g ", acf[i]);
	}
	printf("\n");

	// acvf_opt : Method 0 General Method
	method = 0;
	printf("\n acvf_opt Method 0 General Method \n");
	acvf_opt(inp, N, method, acf, M);

	for (i = 0; i < M; ++i) {
		printf("%g ", acf[i]);
	}

	printf("\n");
	// acvf_opt : Method 1 FFT
	method = 1;
	printf("\n acvf_opt Method 1 FFT \n");
	acvf_opt(inp, N, method, acf, M);

	for (i = 0; i < M; ++i) {
		printf("%g ", acf[i]);
	}
	printf("\n \n");

	// Calculate Autocorrelation from already caluclated Covariance function
	printf("Autocorrelation \n");
	acvf2acf(acf, M);

	for (i = 0; i < M; ++i) {
		printf("%g ", acf[i]);
	}
	printf("\n");

	free(inp);
	free(acf);
	return 0;

}
