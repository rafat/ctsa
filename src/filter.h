// SPDX-License-Identifier: BSD-3-Clause
/*
 * filter.h
 *
 *  Created on: May 26, 2013
 *      Author: Rafat Hussain
 */

#ifndef FILTER_H_
#define FILTER_H_

#include "conv.h"

#ifdef __cplusplus
extern "C" {
#endif

double norm2(double *array, int N);

void conv(double *sig, int len_sig, double *filt, int len_filt, double *oup);

void convx(double *sig, int len_sig, double *filt, int len_filt, double *oup);

void filter(double *sig, int len_sig, double *filt, int len_filt, double *oup);

void mafilter(double *sig,int N, int window,double *oup);

void mafilter2(double *sig,int N, int window,double *oup);

void expfilter(double *sig,int N, double alpha,double *oup);

void mafilter_wt(double *sig,int N, double *weights, int window,double *oup);

#ifdef __cplusplus
}
#endif

#endif /* FILTER_H_ */
