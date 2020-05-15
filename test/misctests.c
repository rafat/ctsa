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

	N = 10;

	x = (double*) malloc(sizeof(double)*N);
	y = (double*) malloc(sizeof(double)*N);

	for(i = 0; i < N; ++i) {
		x[i] = (double) (i+1);
		y[i] = x[i] *x[i];
	}

	Nout = 50;

	xout = (double*) malloc(sizeof(double)*Nout);
	yout = (double*) malloc(sizeof(double)*Nout);

	linspace(xout,Nout,x[0],x[N-1]);

	/* x and y should be sorted in ascending order. x has all unique values
		xout contains equally spaced Nout values from lowest x[0]
		to highest x[N-1].Outputs y[out] are Nout interpolated values.
	*/

	approx(x,y,N,xout,yout,Nout);

	mdisplay(xout,1,Nout);
	mdisplay(yout,1,Nout);

	free(x);
	free(y);
	free(xout);
	free(yout);
}

void arraytest() {
	double x[10] = {0.7,8.2,-21.4,5,8.9,-32.3,7.1,-4.3,1.1,1.2};

	double min,max;

	arrayminmax(x,10,&min,&max);

	printf("MIN %g MAX %g \n",min,max);
}

void interpolatetest() {
	int N;
	double ylo,yhi,val,out;

	double tablei[6] = { -4.38, -4.15, -4.04, -3.99, -3.98, -3.96};
	double tablet[6] = {25, 50, 100, 250, 500, 100000};

	val = 2000;
	N = 6;

	arrayminmax(tablei,N,&ylo,&yhi);


	out = interpolate_linear(tablet,tablei,N,val);

	printf("OUT %g \n",out);


}

void interpolatetest2() {
	int N,i;
	double ylo,yhi,val,out;

	double tablei[4] = { 0.216, 0.176, 0.146, 0.119};
	double tablet[4] = {0.01, 0.025, 0.05, 0.10};
	int pos[4] = {0,0,0,0};

	val = 0.059;
	N = 4;

	arrayminmax(tablei,N,&ylo,&yhi);


	out = interpolate_linear(tablet,tablei,N,val);

	printf("OUT %g \n",out);

	sort1d_ascending(tablei,4,pos);

	for(i = 0; i < 4; ++i) {
		printf(" %d ",pos[i]);
	}



}

void urdf_test() {
	int N,k;
	double stat, pval;

	double x[100] = { -0.5604756, -0.2301775, 1.558708, 0.07050839, 0.1292877, 1.715065, 0.4609162, -1.265061, -0.6868529, -0.445662,
		1.224082, 0.3598138, 0.4007715, 0.1106827, -0.5558411, 1.786913, 0.4978505, -1.966617, 0.7013559, -0.4727914, -1.067824, -0.2179749,
		-1.026004, -0.7288912, -0.6250393, -1.686693, 0.837787, 0.1533731, -1.138137, 1.253815, 0.4264642, -0.2950715, 0.8951257, 0.8781335, 0.8215811,
		0.6886403, 0.5539177, -0.06191171, -0.3059627, -0.380471, -0.694707, -0.2079173, -1.265396, 2.168956, 1.207962, -1.123109, -0.4028848, -0.4666554,
		0.7799651, -0.08336907, 0.2533185, -0.02854676, -0.04287046, 1.368602, -0.225771, 1.516471, -1.548753, 0.5846137, 0.1238542, 0.2159416, 0.3796395,
		-0.5023235, -0.3332074, -1.018575, -1.071791, 0.3035286, 0.4482098, 0.05300423, 0.9222675, 2.050085, -0.4910312, -2.309169, 1.005739, -0.7092008,
		-0.6880086, 1.025571, -0.284773, -1.220718, 0.1813035, -0.1388914, 0.005764186, 0.3852804, -0.37066, 0.6443765, -0.2204866, 0.331782, 1.096839,
		0.4351815, -0.3259316, 1.148808, 0.9935039, 0.548397, 0.2387317, -0.6279061, 1.360652, -0.6002596, 2.187333, 1.532611, -0.2357004, -1.026421 };
	
	N = 100;
	k = -1;
	const char* alternative = "stationary";

	ur_df(x,N,alternative,&k,&stat,&pval);

	printf("Lag %d Stats %g PVAL %g \n",k,stat,pval);
}

void urkpss_test() {
	int N,k,lshort;
	double stat, pval;

	double x[100] = { -0.5604756, -0.2301775, 1.558708, 0.07050839, 0.1292877, 1.715065, 0.4609162, -1.265061, -0.6868529, -0.445662,
		1.224082, 0.3598138, 0.4007715, 0.1106827, -0.5558411, 1.786913, 0.4978505, -1.966617, 0.7013559, -0.4727914, -1.067824, -0.2179749,
		-1.026004, -0.7288912, -0.6250393, -1.686693, 0.837787, 0.1533731, -1.138137, 1.253815, 0.4264642, -0.2950715, 0.8951257, 0.8781335, 0.8215811,
		0.6886403, 0.5539177, -0.06191171, -0.3059627, -0.380471, -0.694707, -0.2079173, -1.265396, 2.168956, 1.207962, -1.123109, -0.4028848, -0.4666554,
		0.7799651, -0.08336907, 0.2533185, -0.02854676, -0.04287046, 1.368602, -0.225771, 1.516471, -1.548753, 0.5846137, 0.1238542, 0.2159416, 0.3796395,
		-0.5023235, -0.3332074, -1.018575, -1.071791, 0.3035286, 0.4482098, 0.05300423, 0.9222675, 2.050085, -0.4910312, -2.309169, 1.005739, -0.7092008,
		-0.6880086, 1.025571, -0.284773, -1.220718, 0.1813035, -0.1388914, 0.005764186, 0.3852804, -0.37066, 0.6443765, -0.2204866, 0.331782, 1.096839,
		0.4351815, -0.3259316, 1.148808, 0.9935039, 0.548397, 0.2387317, -0.6279061, 1.360652, -0.6002596, 2.187333, 1.532611, -0.2357004, -1.026421 };
	
	N = 100;
	k = 0;
	const char* type = "Level";
	lshort = 1;

	ur_kpss(x,N,type,lshort,&k, &stat,&pval);

	printf("Lag %d Stats %g PVAL %g \n",k,stat,pval);
}

void urpp_test() {
	int N,k,lshort;
	double stat, pval;

	double x[100] = { -0.5604756, -0.2301775, 1.558708, 0.07050839, 0.1292877, 1.715065, 0.4609162, -1.265061, -0.6868529, -0.445662,
		1.224082, 0.3598138, 0.4007715, 0.1106827, -0.5558411, 1.786913, 0.4978505, -1.966617, 0.7013559, -0.4727914, -1.067824, -0.2179749,
		-1.026004, -0.7288912, -0.6250393, -1.686693, 0.837787, 0.1533731, -1.138137, 1.253815, 0.4264642, -0.2950715, 0.8951257, 0.8781335, 0.8215811,
		0.6886403, 0.5539177, -0.06191171, -0.3059627, -0.380471, -0.694707, -0.2079173, -1.265396, 2.168956, 1.207962, -1.123109, -0.4028848, -0.4666554,
		0.7799651, -0.08336907, 0.2533185, -0.02854676, -0.04287046, 1.368602, -0.225771, 1.516471, -1.548753, 0.5846137, 0.1238542, 0.2159416, 0.3796395,
		-0.5023235, -0.3332074, -1.018575, -1.071791, 0.3035286, 0.4482098, 0.05300423, 0.9222675, 2.050085, -0.4910312, -2.309169, 1.005739, -0.7092008,
		-0.6880086, 1.025571, -0.284773, -1.220718, 0.1813035, -0.1388914, 0.005764186, 0.3852804, -0.37066, 0.6443765, -0.2204866, 0.331782, 1.096839,
		0.4351815, -0.3259316, 1.148808, 0.9935039, 0.548397, 0.2387317, -0.6279061, 1.360652, -0.6002596, 2.187333, 1.532611, -0.2357004, -1.026421 };
	
	N = 100;
	k = 0;
	const char* type = "Z(t_alpha)";
	const char* alternative = "stationary";
	lshort = 1;

	ur_pp(x,N,alternative,type,lshort,&k, &stat,&pval);

	printf("Lag %d Stats %g PVAL %g \n",k,stat,pval);
}

int main() {
    //errortests();
	//llstest();
	//urtest();
	//approxtest();
	//arraytest();
	//interpolatetest();
	urdf_test();
	//urkpss_test();
	//interpolatetest2();
	urpp_test();
    return 0;
}