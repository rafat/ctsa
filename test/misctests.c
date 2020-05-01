#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/errors.h"
#include "../src/unitroot.h"

void llsptest() {
	int M, N,p,ret,i;
	double *b,*x;
	double val,val2;

	M = 6;
	N = 3;
	p = 2;

	b = (double*)malloc(sizeof(double)* M * p);
	x = (double*)malloc(sizeof(double)* N * p);

	//double A[9] = { 2.8, 1.7, 3, 1, -2, 4, -3, 9, 7 };
	double A[18] = { 100.0,10.0,1.0,104.04,10.2,1.0,108.16,10.4,1.0,112.36,10.6,1.0,116.64,10.8,1.0,121.0,11.0,1.0 };

	val = 0.0;
	val2 = 3.0;
	for (i = 0; i < M; ++i) {
		b[2*i] = val*val / 10;
		b[2*i+1] = 2 * val*val - 1.0;
		val += 0.2;
	}

	mdisplay(b, M, p);
	//ret = lls_qr_p(A, b, M, N,p, x);
	//ret = lls_svd(A, b, M, N, x);

	mdisplay(x, N, p);

	for (i = 0; i < N*p; ++i) {
		x[i] = 0;
	}

	mdisplay(A, M, N);
	mdisplay(b, M, p);
	//ret = lls_svd_p(A, b, M, N,p, x);
	mdisplay(x, N, p);
	printf("RET %d \n", ret);
	
	free(b);
	free(x);
}

void llstest() {
    int M, N,ret,i;
	double *x, *b;

    double val,val2;

	M = 6;
	N = 3;

    b = (double*)malloc(sizeof(double)* M);
	x = (double*)malloc(sizeof(double)* N);

	//double A[9] = { 2.8, 1.7, 3, 1, -2, 4, -3, 9, 7 };
	double A[18] = { 100.0,10.0,1.0,104.04,10.2,1.0,108.16,10.4,1.0,112.36,10.6,1.0,116.64,10.8,1.0,121.0,11.0,1.0 };
	//double A[6] = {1,2,0,0,4,3};

    val = 0.0;
	val2 = 3.0;
	for (i = 0; i < M; ++i) {
		//b[i] = val*val / 10;
		b[i] = 2 * val*val - 1.0;
		val += 0.2;
	}

	printf("B : \n");
	mdisplay(b, M, 1);

	ret = lls_normal(A, b, M, N, x);

	printf("\n X (Normal) : \n");
	mdisplay(x, N, 1);

	for (i = 0; i < N; ++i) {
		x[i] = 0;
	}

	ret = lls_qr(A, b, M, N, x);

	printf("\n X (QR) : \n");
	mdisplay(x, N, 1);

	for (i = 0; i < N; ++i) {
		x[i] = 0;
	}

	ret = lls_svd2(A, b, M, N, x);

	printf("\n X (SVD2) : \n");
	mdisplay(x, N, 1);

	for (i = 0; i < N; ++i) {
		x[i] = 0;
	}

	ret = lls_svd(A, b, M, N, x);

	printf("\n X (SVD) : \n");
	mdisplay(x, N, 1);

	free(b);
	free(x);
	
}

void errortests() {
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
}

void urtest() {
	int N;

	N = 10;

	double y[10] = {1,2,3,4,5,6,7,8,9,10};

	//ur_df(y,N,NULL,2,NULL);
}

void approxtest() {
	int i,N, Nout;
	double *x,*y,*xout,*yout;
	double ylo,yhi;

	N = 10;

	x = (double*) malloc(sizeof(double)*N);
	y = (double*) malloc(sizeof(double)*N);

	for(i = 0; i < N; ++i) {
		x[i] = (double) (i+1);
		y[i] = x[i] *x[i];
	}

	ylo = y[0];
	yhi = y[N-1];

	Nout = 50;

	xout = (double*) malloc(sizeof(double)*Nout);
	yout = (double*) malloc(sizeof(double)*Nout);

	linspace(xout,Nout,x[0],x[N-1]);

	approx(x,y,N,xout,yout,Nout,ylo,yhi);

	mdisplay(xout,1,Nout);
	mdisplay(yout,1,Nout);

	free(x);
	free(y);
	free(xout);
	free(yout);
}

int main() {
    //errortests();
	//llstest();
	//urtest();
	approxtest();
    return 0;
}