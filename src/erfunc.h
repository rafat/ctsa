/*
 * erfunc.h
 *
 *  Created on: Jun 16, 2013
 *      Author: USER
 */

#ifndef ERFUNC_H_
#define ERFUNC_H_


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PIVAL 3.14159265358979323846264338327950288

#define XINFVAL 1.79e+308

#define XNINFVAL 2.2251e-308

#ifndef DBL_MAX_EXP
#define DBL_MAX_EXP   +1024
#endif

#ifndef DBL_MIN_EXP
#define DBL_MIN_EXP   -1021
#endif

double erf__(double x);

double erfc1__(int ind, double x);

double erf(double x);

double erfc(double x);

double erfcx(double x);

double erfinv(double x);

double erfcinv(double x);


#ifdef __cplusplus
}
#endif



#endif /* ERFUNC_H_ */
