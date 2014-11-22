/*
 * neldormead.h
 *
 *  Created on: Jan 5, 2014
 *      Author: HOME
 */

#ifndef NELDERMEAD_H_
#define NELDERMEAD_H_

#include "brent.h"

#ifdef __cplusplus
extern "C" {
#endif

int nel_min(custom_function *funcpt,double *xc,int N,double *dx,double fsval,int MAXITER,int *niter,
		double eps,double *xf);

#ifdef __cplusplus
}
#endif


#endif /* NELDORMEAD_H_ */
