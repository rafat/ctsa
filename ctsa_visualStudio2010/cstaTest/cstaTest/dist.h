/*
 * dist.h
 *
 *  Created on: Jun 11, 2013
 *      Author: USER
 */

#ifndef DIST_H_
#define DIST_H_

#include <stdio.h>
#include <stdlib.h>
#include "erfunc.h"

#ifdef __cplusplus
extern "C" {
#endif


double fix(double x);

double normf(double x);

double eps(double x);

double log1p(double x);

double recfact(double N);

double factorial(int N);

double r8_max ( double x, double y );

double gamma(double x);

double gamma_log(double x);

double beta(double a, double b);

double beta_log(double a, double b);

double pgamma(double x, double a); //Incomplete Gamma (Lower Tail)

double qgamma(double x, double a); ////Incomplete Gamma (Upper Tail)

double bfrac(double x,double a, double b, int* ctr); // Continued Fractions Algorithm (ibeta)

double ibeta_appx(double x, double a , double b); // Incomplete Beta Approximation

double ibeta(double x, double a , double b); //Incomplete Beta

double ibetac(double x, double a , double b); //Incomplete Beta Complement

double ibetad(double x, double a , double b); //Incomplete Beta Derivative

double betapdf(double x, double a , double b);

double betacdf(double x, double a , double b);

double betainv(double alpha, double p, double q); // Beta Inverse

#ifdef __cplusplus
}
#endif



#endif /* DIST_H_ */
