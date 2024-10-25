// SPDX-License-Identifier: BSD-3-Clause
/*
 * pdist.h
 *
 *  Created on: Jun 19, 2013
 *      Author: Rafat Hussain
 */

#ifndef PDIST_H_
#define PDIST_H_

#include "dist.h"

#ifdef __cplusplus
extern "C" {
#endif

// Probability Distributions (pdf,cdf and inverse cdf)

// 1. Normal Distribution

double normalpdf(double x, double mu, double sigma);

double normalcdf(double x, double mu, double sigma);

double normalinv(double p, double mu, double sigma);

// 2. Student's T Distribution

double tpdf(double t, int df);

double tcdf(double t, int df);

double tinv_appx(double p, int df);

double tinv(double p, int df);

// 3. F Distribution

double fpdf(double x, int k1,int k2);

double fcdf(double x, int k1,int k2);

double finv(double p, int k1,int k2);

// 4. Gamma Distribution

double gammapdf(double x, double k, double th);

double gammacdf(double x, double k, double th);

double gammainv(double p, double k, double th);

// 5. Chi-squared Distribution

double chipdf(double x, int df);

double chicdf(double x, int df);

double chiinv(double p, int df);

#ifdef __cplusplus
}
#endif


#endif /* PDIST_H_ */
