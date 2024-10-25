// SPDX-License-Identifier: LGPL-3.0-or-later
// see https://people.math.sc.edu/Burkardt/f_src/brent/brent.html
#ifndef BRENT_H_
#define BRENT_H_

#include "secant.h"


#ifdef __cplusplus
extern "C" {
#endif


double brent_zero(custom_funcuni *funcuni,double a, double b, double tol, double eps);

double brent_local_min(custom_funcuni *funcuni, double a, double b, double t, double eps, double *x);

#ifdef __cplusplus
}
#endif

#endif /* SECANT_H_ */
