/*
 * stats.c
 *
 *  Created on: May 13, 2013
 *      Author: Rafat Hussain
 */

#include "stats.h"

double mean(double* vec, int N) {
	int i;
	double m;
	m = 0.0;

	for (i = 0; i < N; ++i) {
		m+= vec[i];
	}
	m = m / N;
	return m;
}

double var(double* vec, int N) {
	double v,temp,m;
	int i;
	v = 0.0;	
	m = mean(vec,N);
	
	for (i = 0; i < N; ++i) {
		temp = vec[i] - m;
		v+= temp*temp;
	}
	
	v = v / N;
	
	return v;
	
}

auto_fft_object auto_fft_init(int N) {
	auto_fft_object obj = NULL;

		
	obj = (auto_fft_object) malloc (sizeof(struct auto_set));
	 
	//obj->clen = npow2(conv_len);
	//obj->clen = conv_len;
	obj->clen = findnexte(2*N);
	//obj->clen = 125;
	obj->ilen = N;
	
	obj->fobj = fft_real_init(obj->clen,1);
	obj->iobj = fft_real_init(obj->clen,-1);
	
	return obj;
}

void autocovar(double* vec,int N, double* acov,int M) {
	double m,temp1,temp2;
	int i,t;
	m = mean(vec,N);
	
	if ( M > N) {
		M = N-1;
		printf("\n Lag is greater than the length N of the input vector. It is automatically set to length N - 1.\n");
		printf("\n The Output Vector only contains N calculated values.");
	} else if ( M < 0) {
		M = 0;
	}
	
	for(i = 0; i < M;i++) {
		acov[i] = 0.0;
		for (t = 0; t < N-i;t++) {
			temp1 = vec[t] - m;
			temp2 = vec[t+i] - m;
			acov[i]+= temp1*temp2;
		}
		acov[i] = acov[i] / N;
		
	}
	
	
}

void autocorr(double* vec,int N,double* acorr, int M) {
	double var;
	int i;
	if (M > N) {
		M = N - 1;
		printf("\n Lag is greater than the length N of the input vector. It is automatically set to length N - 1.\n");
		printf("\n The Output Vector only contains N calculated values.");
	}
	else if (M < 0) {
		M = 0;
	}
	autocovar(vec,N,acorr,M);
	var = acorr[0];
	acorr[0] = 1.0;
	
	for(i = 1; i < M; i++) {
		acorr[i] = acorr[i]/var;
	}
	
}

void autocovar_fft(auto_fft_object obj,double* vec,double* acov, int M) {
	int i,N,L,N2;
	double m;
	fft_type* a;
	fft_type* b;
	fft_data* ao;
	fft_data* bo;
	
	N = obj->clen;
	L = obj->ilen;
	N2 = N * L;
	if (M > N) {
		M = N - 1;
		printf("\n Lag is greater than the length N of the input vector. It is automatically set to length N - 1.\n");
		printf("\n The Output Vector only contains N calculated values.");
	}
	else if (M < 0) {
		M = 0;
	}
	
	a = (fft_type*) malloc (sizeof(fft_data) * N);
	b = (fft_type*) malloc (sizeof(fft_data) * N);
	ao = (fft_data*) malloc (sizeof(fft_data) * N);
	bo = (fft_data*) malloc (sizeof(fft_data) * N);
	
	m = mean(vec,L);
	
	for (i = 0; i < N;i++) {
		if (i < L) {
			a[i] = vec[i] - m;
		} else {
			a[i] = 0.0;
		}
	
	}
	
	fft_r2c_exec(obj->fobj,a,ao);
	
	for (i = 0; i < N;i++) {
		bo[i].re = ao[i].re * ao[i].re + ao[i].im * ao[i].im;
		bo[i].im = 0.0;
	}
	
	fft_c2r_exec(obj->iobj,bo,b);
	
	for (i = 0; i < M;i++) {
		acov[i] = b[i]/N2;
	}
	
	free(a);
	free(b);
	free(ao);
	free(bo);

}

void autocorr_fft(auto_fft_object obj,double* vec,double* acorr, int M) {
	double var;
	int i,N;
	N = obj->clen;
	if (M > N) {
		M = N - 1;
		printf("\n Lag is greater than the length N of the input vector. It is automatically set to length N - 1.\n");
		printf("\n The Output Vector only contains N calculated values.");
	}
	else if (M < 0) {
		M = 0;
	}
	autocovar_fft(obj,vec,acorr,M);
	var = acorr[0];
	acorr[0] = 1.0;
	
	for(i = 1; i < M; i++) {
		acorr[i] = acorr[i]/var;
	}
	
}

void free_auto(auto_fft_object object) {
	free(object);
}
