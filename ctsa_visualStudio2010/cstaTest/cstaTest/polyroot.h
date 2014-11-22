/*
* polyroot.h
*
*  Created on: August 21, 2014
*      Author: Rafat Hussain
*/

#ifndef POLYROOT_H_
#define POLYROOT_H_

#include <float.h>
#include "stats.h"
#include "matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

int cpoly(double *OPR, double *OPI, int DEGREE, double *ZEROR, double *ZEROI);

int polyroot(double *coeff, int DEGREE, double *ZEROR, double *ZEROI); // Find roots of a real polynomial

int cpolyroot(double *rcoeff,double *icoeff, int DEGREE, double *ZEROR, double *ZEROI);// Find roots of a complex polynomial

#ifdef __cplusplus
}
#endif



#endif /* POLYROOT_H_ */