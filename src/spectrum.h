/*
 * spectrum.h
 *
 *  Created on: May 18, 2013
 *      Author: Rafat Hussain
 */

#ifndef SPECTRUM_H_
#define SPECTRUM_H_

#include "stats.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct spectrum_set* spectrum_fft_object;

spectrum_fft_object spectrum_fft_init(int N,int Nfft);

struct spectrum_set{
	fft_object fobj;
	fft_object fobj2;
	int len;
	int lenfft;
};

int isOdd(int N);

void spectrum_shift(double* spec,int N);

void periodogram(spectrum_fft_object obj,double* vec,double* spec, double *freq, int side);

void psd(spectrum_fft_object obj,double* vec,double* spec, double *freq, int side);

void psd_autocovar(auto_fft_object obj,double *vec,double* spec, double *freq,int Nfft, int side);

void free_spectrum(spectrum_fft_object object);

#ifdef __cplusplus
}
#endif

#endif /* SPECTRUM_H_ */
