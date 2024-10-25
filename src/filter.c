// SPDX-License-Identifier: BSD-3-Clause
/*
 * filter.c
 *
 *  Created on: May 26, 2013
 *      Author: Rafat Hussain
 */

#include "filter.h"

double norm2(double *array, int N) {
	int i;
	double nrm;
	
	nrm = 0.0;
	for (i = 0; i < N; ++i) {
		nrm += array[i]*array[i];
	}
	nrm = sqrt(nrm);
	return nrm;
}

void conv(double *sig, int len_sig, double *filt, int len_filt, double *oup) {
	if (len_sig >= len_filt) {
		conv_direct(sig,len_sig, filt, len_filt,oup);
	} else {
		conv_direct(filt,len_filt, sig, len_sig,oup);
	}
}

void convx(double *sig, int len_sig, double *filt, int len_filt, double *oup) {
	if (len_sig >= len_filt) {
		conv_directx(sig,len_sig, filt, len_filt,oup);
	} else {
		conv_directx(filt,len_filt, sig, len_sig,oup);
	}
}

void filter(double *sig, int len_sig, double *filt, int len_filt, double *oup) {
	int N,i;
	double *temp;
	/*
	 *
	 * filter returns the output which is of same length as the input
	 * unlike conv function
	 *
	 */
	
	N = len_sig + len_filt - 1;
	temp = (double*) malloc (sizeof(double) * N);
	
	conv_direct(sig,len_sig, filt, len_filt,temp);
	
	for (i = 0; i < len_sig; ++i) {
		oup[i] = temp[i];
	}
	
	free(temp);
	
	
}

void mafilter(double *sig,int N, int window,double *oup) {
	int odd,q,N2,i,j;
	double sum;
	double *temp;
	/*
	 * mafilter is a moving average filter that smoothes the signal using
	 * a smoothing window. In this implementation, the signal is extended
	 * by the length of the window. Length floor(window/2) extension is carried
	 * out at the beginning of the signal by setting values of the signal
	 * to sig[0]. The rest of the extension is done at the tail of the signal
	 * by setting it equal to sig[N-1].
	 *
	 */
	
	odd = window - ((window/2) * 2);
	
	if (odd) {
		q = (window - 1 ) / 2; 
	} else {
		q = window / 2;
	}
	
	N2 = N + (2 * q);
	
	temp = (double*) malloc (sizeof(double) * N2);
	
	for (i = 0; i < q; ++i) {
		temp[i] = sig[0];
	}
	
	for (i = q; i < N + q; ++i) {
		temp[i] = sig[i - q];
	}
	
	for (i = N + q; i < N2; ++i) {
		temp[i] = sig[N-1];
	}
	
	if (odd) {
		for (i = q; i < N+q; ++i) {
			sum = 0.0;
			for (j = i - q; j < i + q + 1;j++) {
				sum += temp[j];
			}
			oup[i-q] = sum /window;
		}  
	} else {
		for (i = q; i < N+q; ++i) {
			sum = 0.0;
			for (j = i - q; j < i + q + 1;j++) {
				if ( j == i - q || j == i + q) {
					sum += 0.5 * temp[j];
				} else {
					sum += temp[j];
				}
			}
			oup[i-q] = sum /window;
		}  
	}
	
	
	free(temp);
}

void mafilter2(double *sig,int N, int window,double *oup) {
	int odd,q,N2,i,j;
	double sum;
	/*
	 * mafilter2 is a smoothing filter that smoothes a signal by using
	 * a window. No signal extension is done in this case. Instead,
	 * first floor(window/2) signal terms and last [window - floor(window/2)]
	 * terms are left untouched. The smoothing is performed over the rest
	 * of the terms.
	 *
	 */
	
	odd = window - ((window/2) * 2);
	
	if (odd) {
		q = (window - 1 ) / 2; 
	} else {
		q = window / 2;
	}
	
	if (odd) {
		for (i = 0; i < q; i++) {
			oup[i] = sig[i];
		}
		for (i = q; i < N-q; ++i) {
			sum = 0.0;
			for (j = i - q; j < i + q + 1;j++) {
				sum += sig[j];
			}
			oup[i] = sum /window;
		}  
		for (i = N-q; i < N; i++) {
			oup[i] = sig[i];
		}
	} else {
		for (i = 0; i < q; i++) {
			oup[i] = sig[i];
		}
		for (i = q; i < N-q; ++i) {
			sum = 0.0;
			for (j = i - q; j < i + q + 1;j++) {
				if ( j == i - q || j == i + q) {
					sum += 0.5 * sig[j];
				} else {
					sum += sig[j];
				}
			}
			oup[i] = sum /window;
		}  
		for (i = N-q; i < N; i++) {
			oup[i] = sig[i];
		}
	}
}

void expfilter(double *sig,int N, double alpha,double *oup) {
	int i;
	double beta;
	
	/* alpha accepts values between 0 and 1.
	 * No error check is implemented presently
	 * so make sure that the correct value of 
	 * alpha is entered.  
	 */ 
	beta = 1.0 - alpha;
	
	oup[0] = sig[0];
	
	for (i = 1; i < N; ++i) {
		oup[i] = alpha * sig[i] + beta * oup[i-1];
	}
}

void mafilter_wt(double *sig,int N, double *weights, int window,double *oup) {
	int odd,q,N2,i,j,ct;
	double sum;
	double *temp;
	
	/*  Accepts Normalized weights. 
	 *  If the weight is not normalized then normalize it by
	 *  dividing all the weights by the sum of all the weights
	 *  [1,1,1,1,1] normalizes to [1,1,1,1,1] / 5 = [0.2,0.2,0.2,0.2,0.2]
	 */ 
	
	odd = window - ((window/2) * 2);
	
	if (odd) {
		q = (window - 1 ) / 2; 
		N2 = N + (2 * q);
	} else {
		q = window / 2;
		N2 = N + (2 * q) - 1;
	}
	
	
	temp = (double*) malloc (sizeof(double) * N2);
	
	for (i = 0; i < q; ++i) {
		temp[i] = sig[0];
	}
	
	for (i = q; i < N + q; ++i) {
		temp[i] = sig[i - q];
	}
	
	for (i = N + q; i < N2; ++i) {
		temp[i] = sig[N-1];
	}
	
	if (odd) {
		for (i = q; i < N+q; ++i) {
			sum = 0.0;
			ct = 0;
			for (j = i - q; j < i + q + 1;j++) {
				sum += weights[ct] * temp[j];
				ct++;
			}
			oup[i-q] = sum;
		}  
	} else {
		for (i = q; i < N+q; ++i) {
			sum = 0.0;
			ct = 0;
			for (j = i - q; j < i + q;j++) {
				sum += weights[ct] * temp[j];
				ct++;
			}
			oup[i-q] = sum;
		}  
	}
	
	
	free(temp);
}

