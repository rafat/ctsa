// SPDX-License-Identifier: BSD-3-Clause
#include "stl.h"

/*

STL routines converted to C using f2c and then modified to integrate with
the rest of the code.

c     
c     from netlib/a/stl: no authorship nor copyright claim in the source;
c     presumably by the authors of 
c     
c     R.B. Cleveland, W.S.Cleveland, J.E. McRae, and I. Terpenning,
c     STL: A Seasonal-Trend Decomposition Procedure Based on Loess, 
c     Statistics Research Report, AT&T Bell Laboratories.
c     
c     Converted to double precision by B.D. Ripley 1999.
c     Indented, goto labels renamed, many goto's replaced by `if then {else}'
c     (using Emacs), many more comments;  by M.Maechler 2001-02.
c     

*/

void stl_(double *y, int *n, int *np, int *ns, int *nt, int *nl, int *isdeg, int *itdeg, int *
	ildeg, int *nsjump, int *ntjump, int *nljump, int *ni, int *no, double *rw, double *season, double *trend, 
	double *work)
{
    /* System generated locals */
    int work_dim1, work_offset, i__1;

    /* Local variables */
    int i__, k, newnl, newnp, newns, newnt;
	int userw;

/* Arg */
/* 	n                   : length(y) */
/* 	ns, nt, nl          : spans        for `s', `t' and `l' smoother */
/* 	isdeg, itdeg, ildeg : local degree for `s', `t' and `l' smoother */
/* 	nsjump,ntjump,nljump: ........     for `s', `t' and `l' smoother */
/*       ni, no              : number of inner and outer (robust) iterations */
/* Var */
    /* Parameter adjustments */
    --trend;
    --season;
    --rw;
    --y;
    work_dim1 = *n + 2 * *np;
    work_offset = 1 + work_dim1;
    work -= work_offset;

    /* Function Body */
    userw = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	trend[i__] = 0.0;
/* L1: */
    }
/* the three spans must be at least three and odd: */
    newns = imax(3,*ns);
    newnt = imax(3,*nt);
    newnl = imax(3,*nl);
    if (newns % 2 == 0) {
	++newns;
    }
    if (newnt % 2 == 0) {
	++newnt;
    }
    if (newnl % 2 == 0) {
	++newnl;
    }
/* periodicity at least 2: */
    newnp = imax(2,*np);
    k = 0;
/* --- outer loop -- robustnes iterations */
L100:
    stlstp_(&y[1], n, &newnp, &newns, &newnt, &newnl, isdeg, itdeg, ildeg, 
	    nsjump, ntjump, nljump, ni, &userw, &rw[1], &season[1], &trend[1],
	     &work[work_offset]);
    ++k;
    if (k > *no) {
	goto L10;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[i__ + work_dim1] = trend[i__] + season[i__];
/* L3: */
    }
    stlrwt_(&y[1], n, &work[work_dim1 + 1], &rw[1]);
    userw = 1;
    goto L100;
/* --- end Loop */
L10:
/*     robustness weights when there were no robustness iterations: */
    if (*no <= 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    rw[i__] = 1.0;
/* L15: */
	}
    }
} 

int stless_(double *y, int *n, int *len, int *ideg, int *njump, int *userw, double *rw, double *ys,
	 double *res)
{
    /* System generated locals */
    int i__1, i__2, i__3;
    double d__1;

    /* Local variables */
    int i__, j, k;
    int ok;
    int nsh;
    double delta;
    int nleft, newnj, nright;

/* Arg */
/* Var */
    /* Parameter adjustments */
    --res;
    --ys;
    --rw;
    --y;

    /* Function Body */
    if (*n < 2) {
	ys[1] = y[1];
	return 0;
    }
/* Computing MIN */
    i__1 = *njump, i__2 = *n - 1;
    newnj = imin(i__1,i__2);
    if (*len >= *n) {
	nleft = 1;
	nright = *n;
	i__1 = *n;
	i__2 = newnj;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    d__1 = (double) i__;
	    stlest_(&y[1], n, len, ideg, &d__1, &ys[i__], &nleft, &nright, &
		    res[1], userw, &rw[1], &ok);
	    if (! ok) {
		ys[i__] = y[i__];
	    }
/* L20: */
	}
    } else {
	if (newnj == 1) {
	    nsh = (*len + 1) / 2;
	    nleft = 1;
	    nright = *len;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (i__ > nsh && nright != *n) {
		    ++nleft;
		    ++nright;
		}
		d__1 = (double) i__;
		stlest_(&y[1], n, len, ideg, &d__1, &ys[i__], &nleft, &nright,
			 &res[1], userw, &rw[1], &ok);
		if (! ok) {
		    ys[i__] = y[i__];
		}
/* L30: */
	    }
	} else {
	    nsh = (*len + 1) / 2;
	    i__2 = *n;
	    i__1 = newnj;
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
		if (i__ < nsh) {
		    nleft = 1;
		    nright = *len;
		} else if (i__ >= *n - nsh + 1) {
		    nleft = *n - *len + 1;
		    nright = *n;
		} else {
		    nleft = i__ - nsh + 1;
		    nright = *len + i__ - nsh;
		}
		d__1 = (double) i__;
		stlest_(&y[1], n, len, ideg, &d__1, &ys[i__], &nleft, &nright,
			 &res[1], userw, &rw[1], &ok);
		if (! ok) {
		    ys[i__] = y[i__];
		}
/* L40: */
	    }
	}
    }
    if (newnj != 1) {
	i__1 = *n - newnj;
	i__2 = newnj;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    delta = (ys[i__ + newnj] - ys[i__]) / (double) newnj;
	    i__3 = i__ + newnj - 1;
	    for (j = i__ + 1; j <= i__3; ++j) {
		ys[j] = ys[i__] + delta * (double) (j - i__);
/* L47: */
	    }
/* L45: */
	}
	k = (*n - 1) / newnj * newnj + 1;
	if (k != *n) {
	    d__1 = (double) (*n);
	    stlest_(&y[1], n, len, ideg, &d__1, &ys[*n], &nleft, &nright, &
		    res[1], userw, &rw[1], &ok);
	    if (! ok) {
		ys[*n] = y[*n];
	    }
	    if (k != *n - 1) {
		delta = (ys[*n] - ys[k]) / (double) (*n - k);
		i__2 = *n - 1;
		for (j = k + 1; j <= i__2; ++j) {
		    ys[j] = ys[k] + delta * (double) (j - k);
/* L55: */
		}
	    }
	}
    }
	return 0;
} /* stless_ */

int stlest_(double *y, int *n, int *len, int *ideg, double *xs, double *ys, int *nleft, int *
	nright, double *w, int *userw, double *rw, int *ok)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2;


    /* Local variables */
    double a, b, c__, h__;
    int j;
    double r__, h1, h9, range;

/* Arg */
/* Var */
    /* Parameter adjustments */
    --rw;
    --w;
    --y;

    /* Function Body */
    range = (double) (*n) - 1.;
/* Computing MAX */
    d__1 = *xs - (double) (*nleft), d__2 = (double) (*nright) - *xs;
    h__ = pmax(d__1,d__2);
    if (*len > *n) {
	h__ += (double) ((*len - *n) / 2);
    }
    h9 = h__ * .999;
    h1 = h__ * .001;
    a = 0.0;
    i__1 = *nright;
    for (j = *nleft; j <= i__1; ++j) {
	r__ = (d__1 = (double) j - *xs, fabs(d__1));
	if (r__ <= h9) {
	    if (r__ <= h1) {
		w[j] = 1.0;
	    } else {
/* Computing 3rd power */
		d__2 = r__ / h__;
/* Computing 3rd power */
		d__1 = 1.f - d__2 * (d__2 * d__2);
		w[j] = d__1 * (d__1 * d__1);
	    }
	    if (*userw) {
		w[j] = rw[j] * w[j];
	    }
	    a += w[j];
	} else {
	    w[j] = 0.0;
	}
/* L60: */
    }
    if (a <= 0.0) {
	*ok = 0;
    } else {
	*ok = 1;
	i__1 = *nright;
	for (j = *nleft; j <= i__1; ++j) {
	    w[j] /= a;
/* L69: */
	}
	if (h__ > 0.f && *ideg > 0) {
	    a = 0.f;
	    i__1 = *nright;
	    for (j = *nleft; j <= i__1; ++j) {
		a += w[j] * (double) j;
/* L73: */
	    }
	    b = *xs - a;
	    c__ = 0.0;
	    i__1 = *nright;
	    for (j = *nleft; j <= i__1; ++j) {
/* Computing 2nd power */
		d__1 = (double) j - a;
		c__ += w[j] * (d__1 * d__1);
/* L75: */
	    }
	    if (sqrt(c__) > range * .001) {
		b /= c__;
		i__1 = *nright;
		for (j = *nleft; j <= i__1; ++j) {
		    w[j] *= b * ((double) j - a) + 1.0;
/* L79: */
		}
	    }
	}
	*ys = 0.0;
	i__1 = *nright;
	for (j = *nleft; j <= i__1; ++j) {
	    *ys += w[j] * y[j];
/* L81: */
	}
    }
    return 0;
} /* stlest_ */

int stlfts_(double *x, int *n, int *np, double *trend, double *work)
{
    /* System generated locals */
    int i__1;
	int c__3 = 3;

    /* Parameter adjustments */
    --work;
    --trend;
    --x;

    /* Function Body */
    stlma_(&x[1], n, np, &trend[1]);
    i__1 = *n - *np + 1;
    stlma_(&trend[1], &i__1, np, &work[1]);
    i__1 = *n - (*np << 1) + 2;
    stlma_(&work[1], &i__1, &c__3, &trend[1]);
    return 0;
} /* stlfts_ */

int stlma_(double *x, int *n, int *len, double *ave)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int i__, j, k, m;
    double v, flen;
    int newn;

/* Moving Average (aka "running mean") */
/* ave(i) := mean(x{j}, j = max(1,i-k),..., min(n, i+k)) */
/*           for i = 1,2,..,n */
/* Arg */
/* Var */
    /* Parameter adjustments */
    --ave;
    --x;

    /* Function Body */
    newn = *n - *len + 1;
    flen = (double) (*len);
    v = 0.0;
    i__1 = *len;
    for (i__ = 1; i__ <= i__1; ++i__) {
	v += x[i__];
/* L3: */
    }
    ave[1] = v / flen;
    if (newn > 1) {
	k = *len;
	m = 0;
	i__1 = newn;
	for (j = 2; j <= i__1; ++j) {
	    ++k;
	    ++m;
	    v = v - x[m] + x[k];
	    ave[j] = v / flen;
/* L7: */
	}
    }
    return 0;
} /* stlma_ */

int stlstp_(double *y, int *n, int *np, int *ns, int *nt, int *nl, int *isdeg, int *itdeg, int 
	*ildeg, int *nsjump, int *ntjump, int *nljump, int *ni, int *userw, double *rw, double *season,
	double *trend, double *work)
{
    /* System generated locals */
    int work_dim1, work_offset, i__1, i__2;

    /* Local variables */
    int i__, j;
	int c_false = 0;

/* Arg */
/* Var */
    /* Parameter adjustments */
    --trend;
    --season;
    --rw;
    --y;
    work_dim1 = *n + 2 * *np;
    work_offset = 1 + work_dim1;
    work -= work_offset;

    /* Function Body */
    i__1 = *ni;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[i__ + work_dim1] = y[i__] - trend[i__];
/* L1: */
	}
	stlss_(&work[work_dim1 + 1], n, np, ns, isdeg, nsjump, userw, &rw[1], 
		&work[(work_dim1 << 1) + 1], &work[work_dim1 * 3 + 1], &work[(
		work_dim1 << 2) + 1], &work[work_dim1 * 5 + 1], &season[1]);
	i__2 = *n + (*np << 1);
	stlfts_(&work[(work_dim1 << 1) + 1], &i__2, np, &work[work_dim1 * 3 + 
		1], &work[work_dim1 + 1]);
	stless_(&work[work_dim1 * 3 + 1], n, nl, ildeg, nljump, &c_false, &
		work[(work_dim1 << 2) + 1], &work[work_dim1 + 1], &work[
		work_dim1 * 5 + 1]);
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    season[i__] = work[*np + i__ + (work_dim1 << 1)] - work[i__ + 
		    work_dim1];
/* L3: */
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[i__ + work_dim1] = y[i__] - season[i__];
/* L5: */
	}
	stless_(&work[work_dim1 + 1], n, nt, itdeg, ntjump, userw, &rw[1], &
		trend[1], &work[work_dim1 * 3 + 1]);
/* L80: */
    }
    return 0;
} /* stlstp_ */

int stlrwt_(double *y, int *n, double *fit, double *rw)
{
    /* System generated locals */
    int i__1;
    double d__1, d__2;

    /* Local variables */
    int i__;
    double r__, c1, c9;
    int mid[2];
    double cmad;
	int c__2 = 2;

/* Robustness Weights */
/* 	rw_i := B( |y_i - fit_i| / (6 M) ),   i = 1,2,...,n */
/* 		where B(u) = (1 - u^2)^2  * 1[|u| < 1]   {Tukey's biweight} */
/* 		and   M := median{ |y_i - fit_i| } */
/* Arg */
/* Var */
    /* Parameter adjustments */
    --rw;
    --fit;
    --y;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rw[i__] = (d__1 = y[i__] - fit[i__], fabs(d__1));
/* L7: */
    }
    mid[0] = *n / 2 + 1;
    mid[1] = *n - mid[0] + 1;
    psort_(&rw[1], *n, mid, c__2);
    cmad = (rw[mid[0]] + rw[mid[1]]) * 3.0;
/*     = 6 * MAD */
    c9 = cmad * .999;
    c1 = cmad * .001;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__ = (d__1 = y[i__] - fit[i__], fabs(d__1));
	if (r__ <= c1) {
	    rw[i__] = 1.f;
	} else if (r__ <= c9) {
/* Computing 2nd power */
	    d__2 = r__ / cmad;
/* Computing 2nd power */
	    d__1 = 1.f - d__2 * d__2;
	    rw[i__] = d__1 * d__1;
	} else {
	    rw[i__] = 0.0;
	}
/* L10: */
    }
    return 0;
} /* stlrwt_ */

int stlss_(double *y, int *n, int *np, int *ns, int *isdeg, int *nsjump, int *userw, double *rw, 
	double *season, double *work1, double *work2, double *work3, double *work4)
{
    /* System generated locals */
    int i__1, i__2, i__3;
	int c__1 = 1;

    /* Local variables */
    int i__, j, k, m;
    int ok;
    double xs;
    int nleft, nright;

/* 	called by stlstp() at the beginning of each (inner) iteration */

/* Arg */
/* Var */
    /* Parameter adjustments */
    --work4;
    --work3;
    --work2;
    --work1;
    --rw;
    --y;
    --season;

    /* Function Body */
    if (*np < 1) {
	return 0;
    }
    i__1 = *np;
    for (j = 1; j <= i__1; ++j) {
	k = (*n - j) / *np + 1;
	i__2 = k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work1[i__] = y[(i__ - 1) * *np + j];
/* L10: */
	}
	if (*userw) {
	    i__2 = k;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		work3[i__] = rw[(i__ - 1) * *np + j];
/* L12: */
	    }
	}
	stless_(&work1[1], &k, ns, isdeg, nsjump, userw, &work3[1], &work2[2],
		 &work4[1]);
	xs = 0.;
	nright = imin(*ns,k);
	stlest_(&work1[1], &k, ns, isdeg, &xs, &work2[1], &c__1, &nright, &
		work4[1], userw, &work3[1], &ok);
	if (! ok) {
	    work2[1] = work2[2];
	}
	xs = (double) (k + 1);
/* Computing MAX */
	i__2 = 1, i__3 = k - *ns + 1;
	nleft = imax(i__2,i__3);
	stlest_(&work1[1], &k, ns, isdeg, &xs, &work2[k + 2], &nleft, &k, &
		work4[1], userw, &work3[1], &ok);
	if (! ok) {
	    work2[k + 2] = work2[k + 1];
	}
	i__2 = k + 2;
	for (m = 1; m <= i__2; ++m) {
	    season[(m - 1) * *np + j] = work2[m];
/* L18: */
	}
/* L200: */
    }
    return 0;
} /* stlss_ */


void stlez_(double *y, int *n, int *np, int *ns, int *isdeg, int *itdeg,int *robust, int *no, 
	double *rw, double *season, double *trend, double *work)
{
    /* System generated locals */
    int work_dim1, work_offset, i__1, i__2;
    double d__1;

    /* Local variables */
    int i__, j, ni, nl, nt;
    double difs, dift, mins, mint, maxs, maxt;
    int ildeg;
    double maxds, maxdt;
    int newnp, newns, nljump, nsjump, ntjump;
	int c_false = 0;
	int c_true = 1;

		// Robust

/* Arg */
/* Var */
    /* Parameter adjustments */
    --trend;
    --season;
    --rw;
    --y;
    work_dim1 = *n + 2 * *np;
    work_offset = 1 + work_dim1;
    work -= work_offset;

    /* Function Body */
    ildeg = *itdeg;
    newns = imax(3,*ns);
    if (newns % 2 == 0) {
	++newns;
    }
    newnp = imax(2,*np);
    nt = newnp * 1.5 / (1 - 1.5 / newns) + .5;
    nt = imax(3,nt);
    if (nt % 2 == 0) {
	++nt;
    }
    nl = newnp;
    if (nl % 2 == 0) {
	++nl;
    }
    if (*robust) {
	ni = 1;
    } else {
	ni = 2;
    }
/* Computing MAX */
    i__1 = 1, i__2 = (int) ((double) newns / 10 + .9);
    nsjump = imax(i__1,i__2);
/* Computing MAX */
    i__1 = 1, i__2 = (int) ((double) nt / 10 + .9);
    ntjump = imax(i__1,i__2);
/* Computing MAX */
    i__1 = 1, i__2 = (int) ((double) nl / 10 + .9);
    nljump = imax(i__1,i__2);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	trend[i__] = 0.0;
/* L2: */
    }
    stlstp_(&y[1], n, &newnp, &newns, &nt, &nl, isdeg, itdeg, &ildeg, &nsjump,
	     &ntjump, &nljump, &ni, &c_false, &rw[1], &season[1], &trend[1], &
	    work[work_offset]);
    *no = 0;
    if (*robust) {
	j = 1;
/*        Loop  --- 15 robustness iterations */
L100:
	if (j <= 15) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		work[i__ + work_dim1 * 6] = season[i__];
		work[i__ + work_dim1 * 7] = trend[i__];
		work[i__ + work_dim1] = trend[i__] + season[i__];
/* L35: */
	    }
	    stlrwt_(&y[1], n, &work[work_dim1 + 1], &rw[1]);
	    stlstp_(&y[1], n, &newnp, &newns, &nt, &nl, isdeg, itdeg, &ildeg, 
		    &nsjump, &ntjump, &nljump, &ni, &c_true, &rw[1], &season[
		    1], &trend[1], &work[work_offset]);
	    ++(*no);
	    maxs = work[work_dim1 * 6 + 1];
	    mins = work[work_dim1 * 6 + 1];
	    maxt = work[work_dim1 * 7 + 1];
	    mint = work[work_dim1 * 7 + 1];
	    maxds = (d__1 = work[work_dim1 * 6 + 1] - season[1], fabs(d__1));
	    maxdt = (d__1 = work[work_dim1 * 7 + 1] - trend[1], fabs(d__1));
	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		if (maxs < work[i__ + work_dim1 * 6]) {
		    maxs = work[i__ + work_dim1 * 6];
		}
		if (maxt < work[i__ + work_dim1 * 7]) {
		    maxt = work[i__ + work_dim1 * 7];
		}
		if (mins > work[i__ + work_dim1 * 6]) {
		    mins = work[i__ + work_dim1 * 6];
		}
		if (mint > work[i__ + work_dim1 * 7]) {
		    mint = work[i__ + work_dim1 * 7];
		}
		difs = (d__1 = work[i__ + work_dim1 * 6] - season[i__], fabs(
			d__1));
		dift = (d__1 = work[i__ + work_dim1 * 7] - trend[i__], fabs(
			d__1));
		if (maxds < difs) {
		    maxds = difs;
		}
		if (maxdt < dift) {
		    maxdt = dift;
		}
/* L137: */
	    }
	    if (maxds / (maxs - mins) < .01 && maxdt / (maxt - mint) < .01) 
		    {
		goto L300;
	    }
/* L151: */
	    ++j;
	    goto L100;
	}
/*        end Loop */
L300:
	;
    } else {
/*     	.not. robust */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    rw[i__] = 1.0;
/* L150: */
	}
    }
} 

int psort_(double *a, int n, int *ind, int ni)
{
    int i__, j, k, l, m, p;
    double t;
    int ij, il[16], jl, iu[16], ju;
    int tt;
    int indl[16], indu[16];


/* Partial Sorting ; used for Median (MAD) computation only */

/* Arg */
/* Var */
    /* Parameter adjustments */
    --a;
    --ind;

    /* Function Body */
    if (n < 0 || ni < 0) {
	return 0;
    }
    if (n < 2 || ni == 0) {
	return 0;
    }
    jl = 1;
    ju = ni;
    indl[0] = 1;
    indu[0] = ni;
    i__ = 1;
    j = n;
    m = 1;
/* Outer Loop */
L161:
    if (i__ < j) {
	goto L10;
    }
/*  _Loop_ */
L166:
    --m;
    if (m == 0) {
	return 0;
    }
    i__ = il[m - 1];
    j = iu[m - 1];
    jl = indl[m - 1];
    ju = indu[m - 1];
    if (! (jl <= ju)) {
	goto L166;
    }
/*     while (j - i > 10) */
L173:
    if (! (j - i__ > 10)) {
	goto L174;
    }
L10:
    k = i__;
    ij = (i__ + j) / 2;
    t = a[ij];
    if (a[i__] > t) {
	a[ij] = a[i__];
	a[i__] = t;
	t = a[ij];
    }
    l = j;
    if (a[j] < t) {
	a[ij] = a[j];
	a[j] = t;
	t = a[ij];
	if (a[i__] > t) {
	    a[ij] = a[i__];
	    a[i__] = t;
	    t = a[ij];
	}
    }
L181:
    --l;
    if (a[l] <= t) {
	tt = a[l];
L186:
	++k;
/* L187: */
	if (! (a[k] >= t)) {
	    goto L186;
	}
	if (k > l) {
	    goto L183;
	}
	a[l] = a[k];
	a[k] = tt;
    }
/* L182: */
    goto L181;
L183:
    indl[m - 1] = jl;
    indu[m - 1] = ju;
    p = m;
    ++m;
    if (l - i__ <= j - k) {
	il[p - 1] = k;
	iu[p - 1] = j;
	j = l;
L193:
	if (jl > ju) {
	    goto L166;
	}
	if (ind[ju] > j) {
	    --ju;
	    goto L193;
	}
	indl[p - 1] = ju + 1;
    } else {
	il[p - 1] = i__;
	iu[p - 1] = l;
	i__ = k;
L200:
	if (jl > ju) {
	    goto L166;
	}
	if (ind[jl] < i__) {
	    ++jl;
	    goto L200;
	}
	indu[p - 1] = jl - 1;
    }
    goto L173;
/*     end while */
L174:
    if (i__ != 1) {
	--i__;
L209:
	++i__;
	if (i__ == j) {
	    goto L166;
	}
	t = a[i__ + 1];
	if (a[i__] > t) {
	    k = i__;
/*           repeat */
L216:
	    a[k + 1] = a[k];
	    --k;
	    if (! (t >= a[k])) {
		goto L216;
	    }
/*           until  t >= a(k) */
	    a[k + 1] = t;
	}
	goto L209;
    }
    goto L161;
/* End Outer Loop */
} /* psort_ */
