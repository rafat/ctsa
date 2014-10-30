/*
 * stats.h
 *
 *  Created on: May 13, 2013
 *      Author: Rafat Hussain
 */

#ifndef STATS_H_
#define STATS_H_

#include "conv.h"
#include "pdist.h"
#include "matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct auto_set* auto_fft_object;

auto_fft_object auto_fft_init(int N);

struct auto_set{
	fft_real_object fobj;
	fft_real_object iobj;
	int ilen;// Length if Input
	int clen;// Length of FFT
};

double mean(double* vec, int N);

double var(double* vec, int N);

// N length of the time series
// M Length of the covariance/correlation vector. The maximum lag is M-1 or N-1 , whichever is smaller.

void autocovar(double* vec,int N,double* acov, int M);

void autocorr(double* vec,int N,double* acorr, int M);

void autocovar_fft(auto_fft_object obj,double* vec,double* acov, int M);

void autocorr_fft(auto_fft_object obj,double* vec,double* acorr, int M);

void free_auto(auto_fft_object object);

#ifdef __cplusplus
}
#endif



#endif /* STATS_H_ */
