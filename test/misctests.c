#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/ctsa.h"
#include "../src/errors.h"
#include "../src/autoutils.h"

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
		printf(" %d %g %g",pos[i],tablei[i],tablei[pos[i]]);
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
	const char* alternative = "stationary";

	ur_df(x,N,alternative,NULL,&stat,&pval);

	printf("Stats %g PVAL %g \n",stat,pval);
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
	const char* type = "level";
	lshort = 1;

	ur_kpss(x,N,type,lshort,NULL, &stat,&pval);

	printf("Stats %g PVAL %g \n",stat,pval);
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

	ur_pp(x,N,alternative,type,lshort,NULL, &stat,&pval);

	printf("Stats %g PVAL %g \n",stat,pval);
}

void decomposetest() {
	int N = 32;
	double x[32] = {-50, 175, 149, 214, 247, 237, 225, 329, 729, 809,
       530, 489, 540, 457, 195, 176, 337, 239, 128, 102, 232, 429, 3,
       98, 43, -141, -77, -13, 125, 361, -45, 184};
	int f = 2;
	int ltrend,lseas,lrandom;
	

	decompose(x,N,f,NULL,"multiplicative",NULL,&ltrend,NULL,&lseas,NULL,&lrandom);
}

void lagstests() {
	double *out;
	int rows,cols;
	int N = 32;
	double x[32] = {-50, 175, 149, 214, 247, 237, 225, 329, 729, 809,
       530, 489, 540, 457, 195, 176, 337, 239, 128, 102, 232, 429, 3,
       98, 43, -141, -77, -13, 125, 361, -45, 184};
	int lags = 3;

	out = genLags(x,N,lags,&rows,&cols);

	mdisplay(out,rows,cols);

	free(out);
}

void ocsbtest() {
	int N = 32;
	int differencing_term;
	double stats,crit;
	double x[32] = {-50, 175, 149, 214, 247, 237, 225, 329, 729, 809,
       530, 489, 540, 457, 195, 176, 337, 239, 128, 102, 232, 429, 3,
       98, 43, -141, -77, -13, 125, 361, -45, 184};
	
	int i;
	const char *method = "aic";

	int f = 2;
	int lags = 3;
	int mlags = 3;

	OCSBtest(x,N,f,mlags,method,&stats,&crit);

	differencing_term = stats > crit ? 1 : 0;

	printf("Differencing Term %d \n",differencing_term);

}

void ocsbtest2() {
	int N = 32;
	double x[32] = {-50, 175, 149, 214, 247, 237, 225, 329, 729, 809,
       530, 489, 540, 457, 195, 176, 337, 239, 128, 102, 232, 429, 3,
       98, 43, -141, -77, -13, 125, 361, -45, 184};
	
	reg_object fit;
	double *tval,*pval;
	int i;
	const char *method = "aic";

	int f = 2;
	int lags = 3;
	int mlags = 3;

	fit = fitOCSB(x,N,f,lags,mlags);

	summary(fit);
	confint(fit);
	anova(fit);
	printf("loglik %g aic %g bic %g aicc %g rank %d \n",fit->loglik,fit->aic,fit->bic,fit->aicc,fit->rank);

	tval = (double*)calloc(fit->p,sizeof(double));
	pval = (double*)calloc(fit->p,sizeof(double));

	zerohyp_val(fit,tval,pval);

	for(i = 0; i < fit->p; ++i) {
		tval[i] = (fit->beta+i)->value/(fit->beta+i)->stdErr;
		printf("B%-25d%-20lf \n",i,tval[i]);
	} 
	
	free_reg(fit);
	free(tval);
	free(pval);
}

void airpassengerstest() {
	int N = 144;
	int differencing_term;
	double stats,crit,critsh;
	double season[1] = {0.0};
	double x[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432};
	
	int i;
	const char *method = "aic";
	int Nseas = 1;
	int f = 12;
	int lags = 3;
	int mlags = 3;

	OCSBtest(x,N,f,mlags,method,&stats,&crit);

	differencing_term = stats > crit ? 1 : 0;

	printf("OCSB Differencing Term %d stats %g crit %g \n",differencing_term,stats,crit);

	critsh = 0.64;

	SHtest(x,N,&f,Nseas,season);

	printf("SH Differencing Term %d season %g crit %g \n",differencing_term,*season,critsh);

}

void ausbeertest() {
	int N = 211;
	int differencing_term;
	double stats,crit,critsh;
	double season[1] = {0.0};
	double x[211] = {284., 213., 227., 308.,
                     262., 228., 236., 320.,
                     272., 233., 237., 313.,
                     261., 227., 250., 314.,
                     286., 227., 260., 311.,
                     295., 233., 257., 339.,
                     279., 250., 270., 346.,
                     294., 255., 278., 363.,
                     313., 273., 300., 370.,
                     331., 288., 306., 386.,
                     335., 288., 308., 402.,
                     353., 316., 325., 405.,
                     393., 319., 327., 442.,
                     383., 332., 361., 446.,
                     387., 357., 374., 466.,
                     410., 370., 379., 487.,
                     419., 378., 393., 506.,
                     458., 387., 427., 565.,
                     465., 445., 450., 556.,
                     500., 452., 435., 554.,
                     510., 433., 453., 548.,
                     486., 453., 457., 566.,
                     515., 464., 431., 588.,
                     503., 443., 448., 555.,
                     513., 427., 473., 526.,
                     548., 440., 469., 575.,
                     493., 433., 480., 576.,
                     475., 405., 435., 535.,
                     453., 430., 417., 552.,
                     464., 417., 423., 554.,
                     459., 428., 429., 534.,
                     481., 416., 440., 538.,
                     474., 440., 447., 598.,
                     467., 439., 446., 567.,
                     485., 441., 429., 599.,
                     464., 424., 436., 574.,
                     443., 410., 420., 532.,
                     433., 421., 410., 512.,
                     449., 381., 423., 531.,
                     426., 408., 416., 520.,
                     409., 398., 398., 507.,
                     432., 398., 406., 526.,
                     428., 397., 403., 517.,
                     435., 383., 424., 521.,
                     421., 402., 414., 500.,
                     451., 380., 416., 492.,
                     428., 408., 406., 506.,
                     435., 380., 421., 490.,
                     435., 390., 412., 454.,
                     416., 403., 408., 482.,
                     438., 386., 405., 491.,
                     427., 383., 394., 473.,
                     420., 390., 410.};
	
	int i;
	const char *method = "aicc";
	int Nseas = 1;
	int f = 4;
	int lags = 3;
	int mlags = 3;

	OCSBtest(x,N,f,mlags,method,&stats,&crit);

	differencing_term = stats > crit ? 1 : 0;

	printf("OCSB Differencing Term %d stats %g crit %g \n",differencing_term,stats,crit);

	critsh = 0.64;

	SHtest(x,N,&f,Nseas,season);

	printf("SH Differencing Term %d season %g crit %g \n",differencing_term,*season,critsh);

}

void sunspotstest() {
	FILE *ifp;
	int differencing_term;
	double stats,crit;
	double temp[3000];
	double *x;
    ifp = fopen("../data/sunspots.txt", "r");
	int N;
	int i;
	const char *method = "aic";

	if (!ifp) {
		printf("Cannot Open File");
		exit(100);
	}
	while (!feof(ifp)) {
		fscanf(ifp, "%lf \n", &temp[i]);
		i++;
	}

	N = i;

	printf("N %d \n",N);

	x = (double*) malloc(sizeof(double)*(N-1));

	for(i = 1; i < N;++i) {
		x[i-1] = temp[i];
	}
	N--;

	int f = 12;
	int lags = 3;
	int mlags = 3;

	OCSBtest(x,N,f,mlags,method,&stats,&crit);

	differencing_term = stats > crit ? 1 : 0;

	printf("Differencing Term %d \n",differencing_term);

	fclose(ifp);
	free(x);


}

void psorttest() {
	int N = 10;
	int NI = 2;
	double a[10] = {0.7,-2.3,2.8,9,-4,7.4,-11.6,4.6,16,0};
	int mid[2] = {0,0};

	psort_(a,N,mid,NI);

	for(int i = 0; i < NI; ++i) {
		printf("%d ",mid[i]);
	}
	mdisplay(a,1,N);
}


void stltest() {
	double x[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432
	};

	int N = 144;
	int f = 12;
	const char* s_window_type = "period";
	double *seasonal,*trend,*remainder;
	int *robust = NULL;

	seasonal = (double*)calloc(N,sizeof(double));
	trend = (double*)calloc(N,sizeof(double));
	remainder = (double*)calloc(N,sizeof(double));

	stl(x,N,f,s_window_type,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,seasonal,trend,remainder);

	mdisplay(seasonal,1,N);
	mdisplay(trend,1,N);
	mdisplay(remainder,1,N);
	free(seasonal);
	free(trend);
	free(remainder);
}

void modstltest() {
	double x[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432
	};

	int N = 144;
	int f = 12;
	double *seasonal,*trend,*remainder;
	int s_window = 13;
	double lambda;

	seasonal = (double*)calloc(N,sizeof(double));
	trend = (double*)calloc(N,sizeof(double));
	remainder = (double*)calloc(N,sizeof(double));

	//stl(x,N,f,s_window_type,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,seasonal,trend,remainder);

	modstl(x,N,f,&s_window,NULL,seasonal,trend,remainder);

	mdisplay(seasonal,1,N);
	mdisplay(trend,1,N);
	mdisplay(remainder,1,N);
	free(seasonal);
	free(trend);
	free(remainder);
}

void modstltest2() {
	int N = 211;
	double x[211] = {284., 213., 227., 308.,
                     262., 228., 236., 320.,
                     272., 233., 237., 313.,
                     261., 227., 250., 314.,
                     286., 227., 260., 311.,
                     295., 233., 257., 339.,
                     279., 250., 270., 346.,
                     294., 255., 278., 363.,
                     313., 273., 300., 370.,
                     331., 288., 306., 386.,
                     335., 288., 308., 402.,
                     353., 316., 325., 405.,
                     393., 319., 327., 442.,
                     383., 332., 361., 446.,
                     387., 357., 374., 466.,
                     410., 370., 379., 487.,
                     419., 378., 393., 506.,
                     458., 387., 427., 565.,
                     465., 445., 450., 556.,
                     500., 452., 435., 554.,
                     510., 433., 453., 548.,
                     486., 453., 457., 566.,
                     515., 464., 431., 588.,
                     503., 443., 448., 555.,
                     513., 427., 473., 526.,
                     548., 440., 469., 575.,
                     493., 433., 480., 576.,
                     475., 405., 435., 535.,
                     453., 430., 417., 552.,
                     464., 417., 423., 554.,
                     459., 428., 429., 534.,
                     481., 416., 440., 538.,
                     474., 440., 447., 598.,
                     467., 439., 446., 567.,
                     485., 441., 429., 599.,
                     464., 424., 436., 574.,
                     443., 410., 420., 532.,
                     433., 421., 410., 512.,
                     449., 381., 423., 531.,
                     426., 408., 416., 520.,
                     409., 398., 398., 507.,
                     432., 398., 406., 526.,
                     428., 397., 403., 517.,
                     435., 383., 424., 521.,
                     421., 402., 414., 500.,
                     451., 380., 416., 492.,
                     428., 408., 406., 506.,
                     435., 380., 421., 490.,
                     435., 390., 412., 454.,
                     416., 403., 408., 482.,
                     438., 386., 405., 491.,
                     427., 383., 394., 473.,
                     420., 390., 410.};
	int f = 4;
	double *seasonal,*trend,*remainder;
	int s_window = 13;
	double lambda;

	seasonal = (double*)calloc(N,sizeof(double));
	trend = (double*)calloc(N,sizeof(double));
	remainder = (double*)calloc(N,sizeof(double));

	//stl(x,N,f,s_window_type,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,seasonal,trend,remainder);

	modstl(x,N,f,&s_window,NULL,seasonal,trend,remainder);

	mdisplay(seasonal,1,N);
	mdisplay(trend,1,N);
	mdisplay(remainder,1,N);
	free(seasonal);
	free(trend);
	free(remainder);
}

void mstltest2() {
	double x[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432
	};

	int N = 144,i,j;
	int Nseas = 1;
	int f[1] = {12};
	int iterate = 2;
	double **seasonal,*trend,*remainder;
	int s_window = 13;
	double lambda;

	seasonal = (double**)calloc(Nseas,sizeof(double*));
	trend = (double*)calloc(N,sizeof(double));
	remainder = (double*)calloc(N,sizeof(double));

	for(i = 0; i < Nseas;++i) {
		seasonal[i] = (double*)calloc(N,sizeof(double));
	}

	mstl(x, N,f,&Nseas,&s_window,NULL,&iterate,seasonal,trend,remainder);

	for(i = 0; i < Nseas;++i) {
		printf("\n Seasonal%d :\n",f[i]);
		for(j = 0; j < N;++j) {
			printf("%g ",seasonal[i][j]);
		}
		printf("\n\n");
	}

	mdisplay(trend,1,N);
	mdisplay(remainder,1,N);
	free(seasonal);
	free(trend);
	free(remainder);
}

void mstltest() {
	FILE *ifp;
	double temp[5000];
	double *x;
    ifp = fopen("../data/taylor.txt", "r");
	int N;
	int i,j;

	if (!ifp) {
		printf("Cannot Open File");
		exit(100);
	}
	while (!feof(ifp)) {
		fscanf(ifp, "%lf \n", &temp[i]);
		i++;
	}

	N = i;

	printf("N %d \n",N);

	x = (double*) malloc(sizeof(double)*N);

	for(i = 0; i < N;++i) {
		x[i] = temp[i];
	}
	int Nseas = 2;
	int f[2] = {48,336};
	int iterate = 2;
	double **seasonal,*trend,*remainder;
	int s_window = 13;
	double lambda;

	seasonal = (double**)calloc(Nseas,sizeof(double*));
	trend = (double*)calloc(N,sizeof(double));
	remainder = (double*)calloc(N,sizeof(double));

	for(i = 0; i < Nseas;++i) {
		seasonal[i] = (double*)calloc(N,sizeof(double));
	}

	mstl(x, N,f,&Nseas,&s_window,NULL,&iterate,seasonal,trend,remainder);

	for(i = 0; i < Nseas;++i) {
		printf("\n Seasonal%d :\n",f[i]);
		for(j = 0; j < N;++j) {
			printf("%g ",seasonal[i][j]);
		}
		printf("\n\n");
	}

	mdisplay(trend,1,N);
	mdisplay(remainder,1,N);

	printf("\n %g %g %g %g \n",x[0],x[1],x[N-2],x[N-1]);
	free(seasonal);
	free(trend);
	free(remainder);
	free(x);
	fclose(ifp);
}

static void nulls(double *x, int *N) {
	if (x == NULL) {
		//*x = 4.5;
		*N = 2;
	}
}

void nulltest() {
	double *x = NULL;
	int N;

	nulls(x,&N);

	printf("N %d ",N);
}

void epstest() {
	printf("MACHEPS %g",macheps()/744.44);
}

void boxcoxtest() {
	double x[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432
	};

	int N = 144;
	double *y,*z;
	double lambda;

	y = (double*) malloc(sizeof(double)*N);
	z = (double*) malloc(sizeof(double)*N);

	lambda = boxcox(x,N,NULL,y);
	
	mdisplay(y,1,N);

	inv_boxcox_eval(y,N,lambda,z);

	mdisplay(z,1,N);

	free(y);
	free(z);
}

void supersmoothertest() {
	double *x, *oup;
	double y[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432
	};
	double span,alpha;
	int periodic,i;

	int N = 144;

	x = (double*) malloc(sizeof(double)*N);
	oup = (double*) malloc(sizeof(double)*N);

	for(i = 0; i < N;++i) {
		x[i] = i;
	}

	periodic = 1;
	span = 0.0;
	alpha = -1.0;

	supsmu(x,N,y,NULL,periodic,span,alpha,oup);

	mdisplay(oup,1,N);

	free(x);
	free(oup);
}

void shtest() {
	int N = 32;
	int differencing_term,Nseas;
	double crit;
	double *season;
	double x[32] = {-50, 175, 149, 214, 247, 237, 225, 329, 729, 809,
       530, 489, 540, 457, 195, 176, 337, 239, 128, 102, 232, 429, 3,
       98, 43, -141, -77, -13, 125, 361, -45, 184};
	
	int i;
	int f = 2;
	
	Nseas = 1;
	crit = 0.64;

	season = (double*) malloc(sizeof(double)*Nseas);

	SHtest(x,N,&f,Nseas,season);

	differencing_term = *season > crit ? 1 : 0;

	printf("stat %g Differencing Term %d \n",*season,differencing_term);

	free(season);
}

void nsdiffstest() {
	double x[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432
	};

	int N = 144;
	int f = 12;
	double alpha = 0.05;
	const char *test = "seas";
	int max_D = 1;
	int D;

	//d = ndiffs(x,N,&alpha,test,&max_d);

	D = nsdiffs(x,N,f,&alpha,test,&max_D);

	printf("%d \n",D);
}

void ndiffstest() {
	/*
	double x[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432
	};

	int N = 144;
	*/
	double x[211] = {284., 213., 227., 308.,
                     262., 228., 236., 320.,
                     272., 233., 237., 313.,
                     261., 227., 250., 314.,
                     286., 227., 260., 311.,
                     295., 233., 257., 339.,
                     279., 250., 270., 346.,
                     294., 255., 278., 363.,
                     313., 273., 300., 370.,
                     331., 288., 306., 386.,
                     335., 288., 308., 402.,
                     353., 316., 325., 405.,
                     393., 319., 327., 442.,
                     383., 332., 361., 446.,
                     387., 357., 374., 466.,
                     410., 370., 379., 487.,
                     419., 378., 393., 506.,
                     458., 387., 427., 565.,
                     465., 445., 450., 556.,
                     500., 452., 435., 554.,
                     510., 433., 453., 548.,
                     486., 453., 457., 566.,
                     515., 464., 431., 588.,
                     503., 443., 448., 555.,
                     513., 427., 473., 526.,
                     548., 440., 469., 575.,
                     493., 433., 480., 576.,
                     475., 405., 435., 535.,
                     453., 430., 417., 552.,
                     464., 417., 423., 554.,
                     459., 428., 429., 534.,
                     481., 416., 440., 538.,
                     474., 440., 447., 598.,
                     467., 439., 446., 567.,
                     485., 441., 429., 599.,
                     464., 424., 436., 574.,
                     443., 410., 420., 532.,
                     433., 421., 410., 512.,
                     449., 381., 423., 531.,
                     426., 408., 416., 520.,
                     409., 398., 398., 507.,
                     432., 398., 406., 526.,
                     428., 397., 403., 517.,
                     435., 383., 424., 521.,
                     421., 402., 414., 500.,
                     451., 380., 416., 492.,
                     428., 408., 406., 506.,
                     435., 380., 421., 490.,
                     435., 390., 412., 454.,
                     416., 403., 408., 482.,
                     438., 386., 405., 491.,
                     427., 383., 394., 473.,
                     420., 390., 410.};

	int N = 211;
	double alpha = 0.05;
	const char *test = "adf";
	const char* type = "trend";
	int max_d = 2;
	int d;

	d = ndiffs(x,N,&alpha,test,type,&max_d);

	printf("%d \n",d);
}

void ur_pptest2() {
	double x[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432
	};

	int N = 144;
	const char *type = "Z-tau";
	const char *model = "constant";
	int lshort = 1;
	double stat,pval,teststat;
	double cval[3] = {0,0,0};
	double cprobs[3] = {0,0,0};
	double auxstat[2] = {0,0};
	int laux;

	ur_pp2(x,N,type,model,lshort,NULL,cval,cprobs,auxstat,&laux,&teststat);

	mdisplay(cval,1,3);
	mdisplay(auxstat,1,laux);
	mdisplay(&teststat,1,1);
}

void ur_dftest2() {
	double x[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432
	};

	int N = 144;
	const char *type = "trend";
	const char *selectlags = "aic";
	int lshort = 1;
	double stat,pval;
	double cval[3] = {0,0,0};
	double cprobs[3] = {0,0,0};
	double teststat[3] = {0,0,0};
	int laux;
	int ltstat;
	int lags = 0;
	int cvrows, cvcols;

	//ur_pp2(x,N,type,model,lshort,NULL,cval,cprobs,auxstat,&laux,&teststat);
	ur_df2(x,N,type,&lags,selectlags,cval,&cvrows,&cvcols,cprobs,teststat,&ltstat);

	mdisplay(cval,cvrows,cvcols);
	//mdisplay(auxstat,1,laux);
	mdisplay(teststat,1,ltstat);
}

void regex1() {
	int N, p;
	double alpha;
	double *res2;
	reg_object fit;

	double varcovar[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	alpha = 0.05;

	double tval3[3] = { 0.0, 0.0, 0.0 };
	double pval3[3] = { 0.0, 0.0, 0.0 };
	/*
	TABLE C.4 PER CAPITA PERSONAL CONSUMPTION EXPENDITURE (PPCE) AND PER CAPITA
	PERSONAL DISPOSABLE INCOME (PPDI) IN THE UNITED STATES, 1956�1970,
	IN 1958 DOLLARS, Gujarati D, Basic Econometrics, 4th Ed. McGraw-Hill
	*/

	N = 144;
	p = 0;
	// p = 3 corresponds to number of coefficients (including the intercept) - b0, b1 , b2
	// YY = b0 + b1 * X1 + b2 * X2 + u
	// where u is residual vector of length N
	// contained in XX vector. where XX = [X1,X2] .
	// YY, X1 and X2 are each of length N = 15
	// X1 = [1839,1844,1831,1881,1883,1910,1969,2016,2126,2239,2336,2404,
	//	2487,2535,2595]
	// X2 =[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]



	res2 = (double*)malloc(sizeof(double)* N);

	double XX[30] = { 1839, 1844, 1831, 1881, 1883, 1910, 1969, 2016, 2126, 2239, 2336, 2404,
		2487, 2535, 2595, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
	double YY[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432
	};

	fit = reg_init(N, p);
	setIntercept(fit,0);// Optional as default value for intercept is 1
	setLLSMethod(fit,"qr");
	regress(fit, NULL, YY, res2, varcovar, alpha);// Perform Regression
	// res2 - residuals vector. varcovar - variance-covariance matrix
	// alpha - Used to determine (1-alpha) * 100 % confidence interval
	// alpha = 0.05 for 95% confidence interval
	mdisplay(res2,1,N);
	summary(fit); // summary of regression
	anova(fit); // ANOVA Table
	confint(fit); //Confidence Intervals of Regression Parameters
	zerohyp_val(fit, tval3, pval3); // Obtain Zerohypothesis t-test values

	double inpx[2] = { 2610, 16 };
	double varx[2] = { 0.0, 0.0 };
	double oupx;

	//oupx = pointpredict(fit, inpx, varcovar, varx); // Fit Values. Returns output
	// for a given (p-1) input vector.
	// varx[0] - Variance of Mean Prediction
	// varx[1] - Variance Of Individual Prediction

	printf("Variance(Mean Pred) %lf , Variance(Indiv. Pred) %lf \n", varx[0], varx[1]);

	free_reg(fit);
	free(res2);
}

void seasdummytest() {
	double x[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432};

	int N = 144;

	int rows, cols;
	int f = 12;
	double *oup;

	oup = seasdummy(x,N,f,&rows,&cols);

	mdisplay(oup,rows,cols);

	free(oup);
}

void sdtests() {
	double x[144] = {112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
        115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
        145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
        171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
        196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
        204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
        242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
        284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
        315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
        340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
        360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
        417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432};

	int N = 144;

	int rows, cols;
	int f = 12;

	//SDtest(x,N,f);
}

void arimatest() {
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
	double temp2[1200];

	ifp = fopen("../data/e6.dat", "r");
	i = 0;
	if (!ifp) {
		printf("Cannot Open File");
		exit(100);
	}
	while (!feof(ifp)) {
		fscanf(ifp, "%lf %lf \n", &temp[i],&temp2[i]);
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
}

void sarimatest() {
	int i, N, d, L;
	double *inp;
	int p, q;
	int s, P, D, Q;
	double *phi, *theta;
	double *PHI, *THETA;
	double *xpred, *amse;
	sarima_object obj;
	p = 0;
	d = 1;
	q = 1;
	s = 12;
	P = 0;
	D = 1;
	Q = 1;


	L = 5;

	phi = (double*)malloc(sizeof(double)* p);
	theta = (double*)malloc(sizeof(double)* q);
	PHI = (double*)malloc(sizeof(double)* P);
	THETA = (double*)malloc(sizeof(double)* Q);

	xpred = (double*)malloc(sizeof(double)* L);
	amse = (double*)malloc(sizeof(double)* L);

	FILE *ifp;
	double temp[1200];
	double temp2[1200];

	ifp = fopen("../data/seriesG.txt", "r");
	//ifp = fopen("../data/itdaily.txt", "r");
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


	obj = sarima_init(p, d, q,s,P,D,Q, N);
	sarima_setMethod(obj, 1); // Method 0 ("MLE") is default so this step is unnecessary. The method also accepts values 1 ("CSS") and 2 ("Box-Jenkins")
	sarima_setCSSML(obj,0); // 0 - Only MLE , 1 - CSS + MLE
	//sarima_setOptMethod(obj, 7);// Method 5 ("BFGS") is default . The method also accepts values 0,1,2,3,4,5,6. Check the documentation for details.
	sarima_exec(obj, inp);
	sarima_summary(obj);
	// Predict the next 5 values using the obtained ARIMA model
	sarima_predict(obj, inp, L, xpred, amse);
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
	sarima_free(obj);
	free(inp);
	free(phi);
	free(theta);
	free(PHI);
	free(THETA);
	free(xpred);
	free(amse);
}

void sarimaxtest() {
	int i, N, d, D, L;
	double *inp;
	int p, q, P, Q, s, r;
	double *phi, *theta;
	double *xreg;
	double *xpred, *amse;
	sarimax_object obj;
	int imean = 1;
	p = 0;
	d = 0;
	q = 0;
	s = 0;
	P = 0;
	D = 0;
	Q = 0;
	r = 0;


	L = 0;

	phi = (double*)malloc(sizeof(double)* p);
	theta = (double*)malloc(sizeof(double)* q);

	xpred = (double*)malloc(sizeof(double)* L);
	amse = (double*)malloc(sizeof(double)* L);

	FILE *ifp;
	double temp[1200];
	double temp2[1200];
	double temp3[1200];

	ifp = fopen("../data/e1m.dat", "r");
	i = 0;
	if (!ifp) {
		printf("Cannot Open File");
		exit(100);
	}
	while (!feof(ifp)) {
		fscanf(ifp, "%lf %lf %lf \n", &temp[i],&temp2[i],&temp3[i]);
		i++;
	}
	N = i;

	inp = (double*)malloc(sizeof(double)* N);
	//wmean = mean(temp, N);
	xreg = (double*)malloc(sizeof(double)* N * 2);

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
		xreg[i] = temp2[i];
		xreg[N+i] = temp3[i];
		//printf("%g \n",inp[i]);
	}


	obj = sarimax_init(p, d, q, P, D, Q, s, r ,imean, N);
	sarimax_setMethod(obj, 0); // Method 0 ("CSS-MLE") is default. The method also accepts values 1 ("MLE") and 2 ("CSS")
	//sarimax_setOptMethod(obj, 5);// Method 7 ("BFGS with More Thuente Line search") is default so this step is unnecessary. The method also accepts values 0,1,2,3,4,5,6. Check the documentation for details.
	sarimax_exec(obj, inp,xreg);
	sarimax_summary(obj);
	// Predict the next 5 values using the obtained ARIMA model
	//arima_predict(obj, inp, L, xpred, amse);
	//printf("\n");
	//printf("Predicted Values : ");
	//for (i = 0; i < L; ++i) {
	//	printf("%g ", xpred[i]);
	//}
	//printf("\n");
	//printf("Standard Errors  : ");
	//for (i = 0; i < L; ++i) {
	//	printf("%g ", sqrt(amse[i]));
	//}
	printf("\n");
	sarimax_free(obj);
	free(inp);
	free(phi);
	free(theta);
	free(xpred);
	free(amse);
	free(xreg);
}

void mainverttest() {
	int ret;
	double ma[3] = {1.0, -0.06331,0.1990};
	double ar[3] = {-1.23,0.4};

	mdisplay(ma,1,3);

	invertroot(3,ma);

	mdisplay(ma,1,3);

	ret = archeck(2,ar);
	printf("ARCHECK %d \n",ret);
}

void refittest() {
	int i, N, d, D, L;
	double *inp;
	int p, q, P, Q, s, r;
	int drift,biasadj,method;
	double *xpred, *amse,*xreg,*newxreg;
	sarimax_wrapper_object obj;
	int imean = 1;
    /*
    Make sure all the parameter values are correct and consistent with other values. eg., if xreg is NULL r should be 0
    or if P = D = Q = 0 then make sure that s is also 0. 
     Recheck the values if the program fails to execute.
    */
	p = 0;
	d = 0;
	q = 0;
	s = 0;
	P = 0;
	D = 0;
	Q = 0;
	r = 2;
	int order[3] = {p,d,q};
	int seasonal[4] = {P,D,Q,s};


	L = 0;

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

	drift = 0;
	biasadj = 0;
	method = 0;

	//obj = sarimax_init(p, d, q, P, D, Q, s, r , N);
	obj = sarimax_wrapper(NULL,inp,N,order,seasonal,xreg,r,drift,imean,NULL,biasadj,method);

    /* setMethod()
    Method 0 ("CSS-MLE") is default. The method also accepts values 1 ("MLE") and 2 ("CSS")
    */

	//sarimax_setMethod(obj, 0); 

    /*sarimax_exec(object, input time series, exogenous time series)
        set exogenous to NULL if deadling only with a univariate time series.
    */
	//sarimax_exec(obj, inp,xreg);
	sarimax_wrapper_summary(obj);
	/* sarimax_predict(sarimax_object obj, double *inp, double *xreg, int L,double *newxreg, double *xpred, double *amse)
        inp - Input Time Series
        xreg - Exogenous Time Series
        L - L point prediction
        newxreg - Exogenous Time Series of dimension r * L where r is the number of exogenus time series and L is the length of each
        xpred - L future values
        amse - MSE for L future values
    */

	sarimax_wrapper_predict(obj, inp, xreg, L, newxreg, xpred, amse);
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

	sarimax_wrapper_free(obj);
	free(inp);
	free(xpred);
	free(amse);
    free(xreg);
    free(newxreg);
}

void myarimatest() {
	int i, N, d, D, L;
	double *inp;
	int p, q, P, Q, s, r;
	int drift,biasadj,method;
	double *xpred, *amse,*xreg,*newxreg;
	myarima_object obj;
	int imean = 1;
    /*
    Make sure all the parameter values are correct and consistent with other values. eg., if xreg is NULL r should be 0
    or if P = D = Q = 0 then make sure that s is also 0. 
     Recheck the values if the program fails to execute.
    */
	p = 0;
	d = 1;
	q = 4;
	s = 0;
	P = 0;
	D = 0;
	Q = 0;
	r = 2;
	int order[3] = {p,d,q};
	int seasonal[4] = {P,D,Q,s};
	const char *ic = "aic";
	int trace = 0;
	int constant = 1;
	int approx = 0;
	double offset = 0;
	int rmethod = 0;


	L = 0;

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

	drift = 1;
	biasadj = 0;
	method = 5;

	//obj = sarimax_init(p, d, q, P, D, Q, s, r , N);
	obj = myarima(inp,N,order,seasonal, constant, ic, trace, approx, offset,xreg, r, &rmethod) ;

    /* setMethod()
    Method 0 ("CSS-MLE") is default. The method also accepts values 1 ("MLE") and 2 ("CSS")
    */

	//sarimax_setMethod(obj, 0); 

    /*sarimax_exec(object, input time series, exogenous time series)
        set exogenous to NULL if deadling only with a univariate time series.
    */
	//sarimax_exec(obj, inp,xreg);
	

	myarima_free(obj);
	free(inp);
	free(xpred);
	free(amse);
    free(xreg);
    free(newxreg);
}

void searchtest() {
	int i,j,p_max,q_max,P_max,Q_max,Order_max,d,D,s,stationary;
	myarima_object obj;
	double *inp,*xreg;
	int N,r = 2;
	const char *ic = "aic";
	int approximation = 0;
	double offset = 0;
	int allowdrift,allowmean,method;

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
	N = i;

	inp = (double*)malloc(sizeof(double)* N);
    xreg = (double*)malloc(sizeof(double)* N * r);

	for (i = 0; i < N; ++i) {
		inp[i] = temp[i];
        xreg[i] = temp1[i];
		xreg[N+i] = temp2[i];
	}
	d = 1;
	D = 0;
	s = 0;
	P_max = 0;
	Q_max = 0;
	p_max = 5;
	q_max = 5;
	Order_max = 5;
	stationary = 0;
	allowdrift = 1;
	allowmean = 1;
	method = 0;

	obj = search_arima(inp,N,d,D,p_max,q_max,P_max,Q_max,Order_max,stationary,s,ic,approximation, xreg,r,offset,allowdrift,allowmean,method);

	printf("ic %g \n",obj->ic);

	free(inp);
	free(xreg);
	myarima_free(obj);
}

int main() {
    //errortests();
	//llstest();
	//urtest();
	//approxtest();
	//arraytest();
	//interpolatetest();
	//urdf_test();
	//urkpss_test();
	//interpolatetest2();
	//urpp_test();
	//decomposetest();
	//lagstests();
	//ocsbtest();
	//ocsbtest2();
	//airpassengerstest();
	//ausbeertest();
	//sunspotstest();
	//psorttest();
	//stltest();
	//nulltest();
	//epstest();
	//boxcoxtest();
	//supersmoothertest();
	//modstltest();
	//modstltest2();
	//mstltest();
	//shtest();
	//ndiffstest();
	//nsdiffstest();
	//ur_pptest2() ;
	//ur_dftest2();
	//regex1();
	//seasdummytest();
	//sdtests();
	//arimatest();
	//sarimatest();
	//sarimaxtest();
	//arimatest();
	//mainverttest();
	//refittest();
	myarimatest();
	searchtest();
    return 0;
}
