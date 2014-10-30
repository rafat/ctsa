/*
============================================================================
Name        : tseries.c
Author      : Rafat Hussain
Version     :
Copyright   :
Description : Hello World in C, Ansi-style
============================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../header/ctsa.h"

int main(void) {
	int i, N, M;
	double *inp, *par;
	int method;

	FILE *ifp;
	double temp[1200];

	ifp = fopen("seriesC.txt", "r");
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
	par = (double*)malloc(sizeof(double)* M);

	// Default Method is Yule-Walker

	printf("\n Default Method : pacf \n");

	pacf(inp, N, par, M);

	for (i = 0; i < M; ++i) {
		printf("%g ", par[i]);
	}
	printf("\n");

	// pacf_opt : Method 0 Yule Walker 
	method = 0;
	printf("\n pacf_opt Method 0 Yule Walker \n");
	pacf_opt(inp, N, method, par, M);

	for (i = 0; i < M; ++i) {
		printf("%g ", par[i]);
	}

	printf("\n");
	// pacf_opt : Method 1 Burg 
	method = 1;
	printf("\n pacf_opt Method 1 Burg \n");
	pacf_opt(inp, N, method, par, M);

	for (i = 0; i < M; ++i) {
		printf("%g ", par[i]);
	}
	printf("\n");

	// pacf_opt : Method 2 MLE (Box-Jenkins) 
	method = 2;
	printf("\n pacf_opt Method 2 MLE (Box-Jenkins) \n");
	pacf_opt(inp, N, method, par, M);

	for (i = 0; i < M; ++i) {
		printf("%g ", par[i]);
	}
	printf("\n");

	free(inp);
	free(par);
	return 0;

}
