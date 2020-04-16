#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/errors.h"

int main() {
    double *predicted, *actual;
    int N,i,M;
    FILE *ifp;
	double temp[1200];

    ifp = fopen("../data/eusm.txt", "r");
	i = 0;
	if (!ifp) {
		printf("Cannot Open File");
		exit(100);
	}
	while (!feof(ifp)) {
		fscanf(ifp, "%lf \n", &temp[i]);
		i++;
	}

    M = i;
    N = 100;


    actual = (double*) malloc(sizeof(double)*N);
    predicted = (double*) malloc(sizeof(double)*N);

    for(i = 0; i < N;++i) {
        predicted[i] = 1716.18;
        actual[i] = temp[200+i];
    }

    printf("Mean Error %g \n",me(predicted,actual,N));
    printf("Mean Squared Error %g \n",mse(predicted,actual,N));
    printf("Root Mean Squared Error %g \n",rmse(predicted,actual,N));
    printf("Mean Absolute Error %g \n",mae(predicted,actual,N));
    printf("Mean Percentage Error %g \n",mpe(predicted,actual,N));
    printf("Mean Absolute Percentage Error %g \n",mape(predicted,actual,N));
    printf("Mean Absolute Scaled Error %g \n",mase(predicted,actual,N,temp,M-N));

    free(actual);
    free(predicted);
    return 0;
}