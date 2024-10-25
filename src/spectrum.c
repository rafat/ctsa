// SPDX-License-Identifier: BSD-3-Clause
/*
 * spectrum.c
 *
 *  Created on: May 18, 2013
 *      Author: Rafat Hussain
 */

#include <stdio.h>
#include "spectrum.h"

spectrum_fft_object spectrum_fft_init(int N,int Nfft) {
	spectrum_fft_object obj = NULL;
	
	obj = (spectrum_fft_object) malloc (sizeof(struct spectrum_set));
	obj->len = N;
	obj->lenfft = Nfft;
	obj->fobj = fft_init(N,1);
	if ( N == Nfft) {
		obj->fobj2 = obj->fobj;
	} else {
		obj->fobj2 = fft_init(Nfft,1);
	}
	
	return obj;
}

int isOdd(int N) {
	return N - ((N/2) * 2);
}

void spectrum_shift(double* spec,int N) {
	int N2,i;
	double temp,temp1;
	
	if (isOdd(N)) {
		N2 = N/2;
		temp1 = spec[0];
		for (i = 0; i < N2; i++) {
			temp = spec[i+1];
			spec[i] = spec[i + N2+1];
			spec[i + N2+1] = temp;
		}
		spec[N2] = temp1;
	
	} else {
		N2 = N/2 - 1;
		temp = spec[N2 + 1];
		temp1 = spec[0]; 
		
		for (i = 0; i < N2; i++) {
			spec[i] = spec[i + N2+2];
			spec[i + N2+1] = spec[i+1];
		}
		spec[N2] = temp1;
		spec[N-1] = temp;
	}
}

void periodogram(spectrum_fft_object obj,double* vec,double* spec, double *freq, int side) {
	int i,N,t,N2;
	double m;
	fft_data *inp;
	fft_data *oup;
	
	N = obj->len;
	
	inp = (fft_data*) malloc (sizeof(fft_data) * N);
	oup = (fft_data*) malloc (sizeof(fft_data) * N);
	
	m = mean(vec,N);
	for (i = 0; i < N;++i) {
		inp[i].re = vec[i];
		inp[i].im = 0.0;
	}
	
	fft_exec(obj->fobj,inp,oup);
	
	if (isOdd(N)) {
		spec[0] = oup[0].re * oup[0].re / N;
		freq[0] = 0.0;
		N2 = N/2;
		for (i = 1; i < N2+1 ;++i) {
			spec[i] = (oup[i].re * oup[i].re + oup[i].im * oup[i].im) / N;
			freq[i] = PI2 * i / N;
		}
		if (side == 2) { 
			for (i = 0; i < N2; i++) {
				spec[N2 + i + 1] = spec[N2 - i];
				freq[N2 + i + 1] = - freq[N2 - i];
			}
		} else {
			for (i = 1; i < N2+1 ;++i) {
				spec[i]*= 2.0;
			}
		}
		
		
	} else {
		spec[0] = oup[0].re * oup[0].re / N;
		freq[0] = 0.0;
		N2 = N/2;
		for (i = 1; i < N2;++i) {
			spec[i] = (oup[i].re * oup[i].re + oup[i].im * oup[i].im) / N;
			freq[i] = PI2 * i / N;
		}
		
		spec[N2] = 0.0;
		freq[N2] = PI2 * N2 / N;
	
		for (t = 0; t < N;t+=2) {
			spec[N2] -= vec[t];
		}
		
		for (t = 1; t < N;t+=2) {
			spec[N2] += vec[t];
		}
		
		spec[N2] = spec[N2] * spec[N2] / N;
		
		if (side == 2) {
			for (i = 1; i < N2; i++) {
				spec[N2 + i] = spec[N2 - i];
				freq[N2 + i] = - freq[N2 - i];
			}
		} else {
			for (i = 1; i < N2;++i) {
				spec[i] *= 2.0;
			}
		}
	}	

}

void psd(spectrum_fft_object obj,double* vec,double* spec, double *freq, int side) {
	int i,N,N2,NO;
	double m,NP;
	fft_data *inp;
	fft_data *oup;
	
	N = obj->lenfft;
	NO = obj->len;
	NP = NO * PI2;
	
	inp = (fft_data*) malloc (sizeof(fft_data) * N);
	oup = (fft_data*) malloc (sizeof(fft_data) * N);
	
	m = mean(vec,N);
	for (i = 0; i < NO;++i) {
		inp[i].re = vec[i];
		inp[i].im = 0.0;
	}
	
	for (i = NO; i < N;++i) {
		inp[i].re = 0.0;
		inp[i].im = 0.0;
	}
	
	fft_exec(obj->fobj2,inp,oup);
	
	if (isOdd(N)) {
		spec[0] = oup[0].re * oup[0].re / NP;
		freq[0] = 0.0;
		N2 = N/2;
		for (i = 1; i < N2+1 ;++i) {
			spec[i] = (oup[i].re * oup[i].re + oup[i].im * oup[i].im) / NP;
			freq[i] = PI2 * i / N;
		}
		if (side == 2) { 
			for (i = 0; i < N2; i++) {
				spec[N2 + i + 1] = spec[N2 - i];
				freq[N2 + i + 1] = - freq[N2 - i];
			}
		} else {
			for (i = 1; i < N2+1 ;++i) {
				spec[i]*= 2.0;
			}
		}
		
		
	} else {
		spec[0] = oup[0].re * oup[0].re / NP;
		freq[0] = 0.0;
		N2 = N/2;
		for (i = 1; i < N2+1;++i) {
			spec[i] = (oup[i].re * oup[i].re + oup[i].im * oup[i].im) / NP;
			freq[i] = PI2 * i / N;
		}
		
		/*
		spec[N2] = 0.0;
		freq[N2] = PI2 * N2 / N;
	
		for (t = 0; t < N;t+=2) {
			spec[N2] -= vec[t];
		}
		
		for (t = 1; t < N;t+=2) {
			spec[N2] += vec[t];
		}
		
		spec[N2] = spec[N2] * spec[N2] / N;
		*/
		if (side == 2) {
			for (i = 1; i < N2; i++) {
				spec[N2 + i] = spec[N2 - i];
				freq[N2 + i] = - freq[N2 - i];
			}
		} else {
			for (i = 1; i < N2;++i) {
				spec[i] *= 2.0;
			}
		}
	}	

}

void psd_autocovar(auto_fft_object obj,double *vec,double* spec, double *freq,int Nfft, int side) {
	int N,f,N2;
	double NP;
	double *acov;
	fft_data *inp;
	fft_data *oup;
	
	N = obj->ilen;
	NP = PI2 * N;
	acov = (double*) malloc (sizeof(double) * N);
	inp = (fft_data*) malloc (sizeof(fft_data) * Nfft);
	oup = (fft_data*) malloc (sizeof(fft_data) * Nfft);
	autocovar_fft(obj,vec,acov,N);
	
	N2 = Nfft/2 + 1;
	fft_object fobj = fft_init(Nfft,1);
	
	for (f = 0; f < N;f++) {
		inp[f].re = acov[f];
		inp[f].im = 0.0;
	} 
	
	for (f = N; f < Nfft;f++) {
		inp[f].re = 0.0;
		inp[f].im = 0.0;
	} 
	
	fft_exec(fobj,inp,oup);
	spec[0] = oup[0].re * oup[0].re / NP;
	freq[0] = 0.0;
	
	for (f = 1; f < N2;++f) {
		spec[f] = (oup[f].re * oup[f].re + oup[f].im * oup[f].im) / NP;
		freq[f] = PI2 * f / Nfft;
	}
	/*
	for (f = 0; f < N2;f++) {
		spec[f] = 0.0;
		for (k = 1; k < N;k++) {
			spec[f]+= acov[k] * cos(PI2 * f * k / Nfft);
		}
		spec[f] *= 2.0;
		spec[f] += acov[0];
		//spec[f] /= PI2;
		freq[f] = PI2 * f / Nfft;
	}
	*/
	
	free_fft(fobj);
	free(acov);
}

void free_spectrum(spectrum_fft_object object) {
	free(object);
}


