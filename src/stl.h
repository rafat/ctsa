// SPDX-License-Identifier: BSD-3-Clause
#ifndef STL_H_
#define STL_H_

#include "polyroot.h"

#ifdef __cplusplus
extern "C" {
#endif

void stl_(double *y, int *n, int *np, int *ns, int *nt, int *nl, int *isdeg, int *itdeg, int *
	ildeg, int *nsjump, int *ntjump, int *nljump, int *ni, int *no, double *rw, double *season, double *trend, 
	double *work);

int stless_(double *y, int *n, int *len, int *ideg, int *njump, int *userw, double *rw, double *ys,
	 double *res);

int stlest_(double *y, int *n, int *len, int *ideg, double *xs, double *ys, int *nleft, int *
	nright, double *w, int *userw, double *rw, int *ok);

int stlfts_(double *x, int *n, int *np, double *trend, double *work);

int stlma_(double *x, int *n, int *len, double *ave);

int stlstp_(double *y, int *n, int *np, int *ns, int *nt, int *nl, int *isdeg, int *itdeg, int 
	*ildeg, int *nsjump, int *ntjump, int *nljump, int *ni, int *userw, double *rw, double *season,
	double *trend, double *work);

int stlrwt_(double *y, int *n, double *fit, double *rw);

int stlss_(double *y, int *n, int *np, int *ns, int *isdeg, int *nsjump, int *userw, double *rw, 
	double *season, double *work1, double *work2, double *work3, double *work4);

void stlez_(double *y, int *n, int *np, int *ns, int *isdeg, int *itdeg,int *robust, int *no, 
	double *rw, double *season, double *trend, double *work);

int psort_(double *a, int n, int *ind, int ni);

#ifdef __cplusplus
}
#endif



#endif /* STL_H_ */