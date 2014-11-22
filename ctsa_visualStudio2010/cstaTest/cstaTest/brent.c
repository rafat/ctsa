#include "brent.h"

/*
 * brent.c
 *
 *  Copyright (c) 2014, Rafat Hussain
 *	License : BSD 3-Clause
 *	See COPYRIGHT for more details
 */

/*
 * Algorithms are C translation of Richard Brent's Algol 60 procedure as published in
 * Algoritms For Minimization Without Derivatives, Richard P. Brent
 * Chapters 3-5
 */ 
double brent_zero(custom_funcuni *funcuni, double ai, double bi, double tol, double eps) {
	double bz;
	double a,b,c,d,e,fa,fb,fc;
	double m,s,etol,fd,p,q,r;
	int iter;
	
	fd = eps;
	
	a = ai;
	b = bi;
	c = a;
	fa = FUNCUNI_EVAL(funcuni,a);
	fb = FUNCUNI_EVAL(funcuni,b);
	fc = fa;
	e = b - a;
	d = e;
	
	if (fabs(fc) < fabs(fb)) {
		a = b;
		b = c;
		c = a;
		fa = fb;
		fb = fc;
		fc = fa;
	}
		
	etol = 2 * fd * fabs(b) + tol;
	m = 0.5 * (c - b);
	
	bz = b;
	iter = 0;
	while (fabs ( m ) >= etol && fb != 0.0) {

		iter++;
		
		if (fabs(e) < etol || fabs(fa) <= fabs(fb)) {
			e = m;
			d = e;
		} else {
			
			s = fb/fa;
			if (a == c) {
				// linear interpolation
				p = 2.0 * m * s;
				q = 1.0 - s;
			} else {
				// inverse quadratic interpolation
				q = fa/fc;
				r = fb/fc;
				p = s * (2 * m * q * (q - r) - (b - a) * (r - 1.0));
				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
			}
			
			if (p > 0.0) {
				q = - q;
			} else {
				p = - p;
			}

			s = e;
			e = d;
			
			if ( 2.0 * p < 3.0 * m * q - fabs( etol * q ) && p < fabs( 0.5 * s * q ) ) {
				d = p /q;
			} else {
				e = m;
				d = e;
			}
			
		}
		a = b;
		fa = fb;
		
		if ( fabs(d) > etol ) {
			b += d;
		} else if ( 0.0 < m ) {
			b += etol;
		} else {
			b -= etol;
		}
		
		fb = FUNCUNI_EVAL(funcuni,b);
		
		bz = b;
		if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) ) {
			c = a;
			fc = fa;
			e = b - a;
			d = e;
		}
      
		// To the end of the loop
		
		if (fabs(fc) < fabs(fb)) {
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		
		etol = 2 * fd * fabs(b) + tol;
		m = 0.5 * (c - b);
		
	}
	
	return bz;
}

double brent_local_min(custom_funcuni *funcuni, double a, double b, double t, double eps, double *x) {
	double c,d,e,m,p,q,r,tol,t2;
	double u,v,w,fu,fv,fw,fx;
	double fd;
	
	fd = eps;
	
	c = (3.0 - sqrt(5.0)) / 2.0;
	*x = a + c * (b - a);
	w = *x; v = w;
	e = 0;
	fx = FUNCUNI_EVAL(funcuni,*x);
	fw = fx; fv = fw;
	
	m = 0.5 * (a + b);
	tol = fd * fabs(*x) + t;
	t2 = 2.0 * tol;
	
	while (fabs(*x - m) > t2 - 0.5 * (b - a)) {
		p = 0; q = 0; r = 0;
		
		if (fabs(e) > tol) {
			r = (*x - w) * (fx - fv);
			q = (*x - v) * (fx - fw);
			p = (*x - v) * q - (*x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0) {
				p = -p;
			} else {
				q = -q;
			}
			r = e;
			e = d;
		}
		
		if (fabs(p) < fabs(0.5 * q * r) && p < q * (a - *x) && p < q * (b - *x)) {
			d = p / q;
			u = *x + d;
			
			if ( (u - a) < t2 || (b - u) < t2) {
				if (*x < m) {
					d = tol;
				} else {
					d = -tol;
				}
			}
		} else {
			if (*x < m) {
				e = b - *x;
			} else {
				e = a - *x;
			}
			d = c * e;
		}
		
		if (fabs(d) >= tol) {
			u = *x + d;
		} else if (d > 0) {
			u = *x + tol;
		} else {
			u = *x - tol;
		}
		
		fu = FUNCUNI_EVAL(funcuni,u);
		
		if (fu <= fx) {
			if (u < *x) {
				b = *x;
			} else {
				a = *x;
			}
			v = w;
			fv = fw;
			w = *x;
			fw = fx;
			*x = u;
			fx = fu;
		} else {
			if (u < *x) {
				a = u;
			} else {
				b = u;
			}
			if (fu <= fw || w == *x) {
				v = w;
				fv = fw;
				w = u;
				fw = fu;
			} else if (fu <= fv || v == *x || v == w) {
				v = u;
				fv = fu;
			}
		}
		
		//End of the loop
		m = 0.5 * (a + b);
		tol = fd * fabs(*x) + t;
		t2 = 2.0 * tol;
		
	}
	
	return fx;
}
