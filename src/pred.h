// SPDX-License-Identifier: BSD-3-Clause
/*
 * pred.h
 *
 *  Created on: Jul 14, 2014
 *      Author: Rafat Hussaint
 */

#ifndef PRED_H_
#define PRED_H_

#include "initest.h"

#ifdef __cplusplus
extern "C" {
#endif

void predictarima(double *zt,int lenzt,int p,int d, int q,double *phi, double *theta,
		double constant,int forlength,double *oup);

#ifdef __cplusplus
}
#endif

#endif /* PRED_H_ */
