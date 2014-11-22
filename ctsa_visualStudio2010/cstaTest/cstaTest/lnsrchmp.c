#include "lnsrchmp.h"

/*
 * lnsrchmp.c
 *
 *  Copyright (c) 2014, Rafat Hussain
 *	License : BSD 3-Clause
 *	See COPYRIGHT for more details
 */

/*
This code is a C translation of Fortran Routine developed by
Argonne National Laboratory. MINPACK Project. June 1983
     Jorge J. More', David J. Thuente

A Copy of the original Fortran routine is available at
http://www.cs.umd.edu/~oleary/LBFGS/FORTRAN/linesearch.f

*/



int grad_fd(custom_function *funcpt, custom_gradient *funcgrad, double *x, int N, double *dx,
		double eps2, double *f) {
	int retval;
	retval = 0;
	if (funcgrad == NULL) {
		//printf("FD Gradient \n");
		retval = grad_calc(funcpt,x,N,dx,eps2,f);
	} else {
		//printf("Analytic gradient \n");
		FUNCGRAD_EVAL(funcgrad,x,N,f);
	}

	return retval;

}

int grad_cd(custom_function *funcpt, custom_gradient *funcgrad, double *x, int N, double *dx,
		double eps3, double *f) {
	int retval;
	retval = 0;
	if (funcgrad == NULL) {
		//printf("FD Gradient \n");
		retval = grad_calc2(funcpt,x,N,dx,eps3,f);
	} else {
		//printf("Analytic gradient \n");
		FUNCGRAD_EVAL(funcgrad,x,N,f);
	}
	return retval;

}

int grad_calc2(custom_function *funcpt, double *x, int N, double *dx, double eps3, double *f) {
	int j,retval;
	double stepsize,stepmax,temp;
	double fp,fm;

	retval = 0;

	for (j = 0; j < N;++j) {
		if (fabs(x[j]) >= 1.0 / fabs(dx[j])) {
			stepmax = x[j];
		}
		else {
			stepmax = signx(x[j]) * 1.0 / fabs(dx[j]);
		}

		stepsize = stepmax * eps3;
		temp = x[j];
		x[j] += stepsize;
		stepsize = x[j] - temp;
		fp = FUNCPT_EVAL(funcpt,x,N);
		if (fp >= DBL_MAX || fp <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			return 15;
		}
		if (fp != fp) {
			printf("Program Exiting as the function returns NaN");
			return 15;
		}
		x[j] = temp - stepsize;
		fm = FUNCPT_EVAL(funcpt,x,N);
		if (fm >= DBL_MAX || fm <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			return 15;
		}
		if (fm != fm) {
			printf("Program Exiting as the function returns NaN");
			return 15;
		}
		f[j] = (fp - fm)/ (2 * stepsize);
		x[j] = temp;
	}

	return retval;
}

int grad_calc(custom_function *funcpt, double *x, int N, double *dx, double eps2, double *f) {
	int i, j,retval;
	double step, fd, stepmax;
	double *xi;

	fd = eps2; // square root of macheps
	retval = 0;
	xi = (double*)malloc(sizeof(double)*N);

	for (i = 0; i < N; ++i) {
		if (fabs(x[i]) >= 1.0 / fabs(dx[i])) {
			stepmax = x[i];
		}
		else {
			stepmax = signx(x[i]) * 1.0 / fabs(dx[i]);
		}
		step = fd * stepmax;
		for (j = 0; j < N; ++j) {
			xi[j] = x[j];
		}
		xi[i] += step;
		f[i] = (FUNCPT_EVAL(funcpt, xi, N) - FUNCPT_EVAL(funcpt, x, N)) / step;
		if (f[i] >= DBL_MAX || f[i] <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			free(xi);
			return 15;
		}
		if (f[i] != f[i]) {
			printf("Program Exiting as the function returns NaN");
			free(xi);
			return 15;
		}
		//xi[i] -= step;
	}

	free(xi);
	return retval;

}

int stopcheck3_mt(double *xi,double *xf,double fx, int N, double fo, double *jac, double *dx, double eps,
		double stoptol, double functol, int retval) {
	int rcode,i;
	double nrm,nrmnx,relfit,num,den,stop0;
	double *scheck;
	rcode = 0;

	scheck = (double*)malloc(sizeof(double)*N);

	if (retval == 3) {
		rcode = 4;
		return rcode;
	}
	if (retval == 15) {
		rcode = 15;
		return rcode;
	}

	nrm = l2norm(jac, N);
	nrmnx = nrm / (double) N;

	if (fabs(fo) < eps) {
		relfit = fabs(fx - fo);
	}
	else {
		relfit = fabs((fx - fo)/fo);
	}

	if (nrmnx < stoptol) {
		rcode = 1; // Successful Convergence
	} else if (relfit < functol) {
		rcode = 6; // Relative fit less than function tolerance
	} else {
		for (i = 0; i < N; ++i) {
			den = 1.0+fabs(xf[i]);
			num = fabs(xf[i] - xi[i]);
			scheck[i] = num / den;
		}
		stop0 = array_max_abs(scheck, N);
		if (stop0 <= stoptol) {
			rcode = 2;
		}
	}

	free(scheck);
	return rcode;
}

int stopcheck2_mt(double fx, int N, double fo, double *jac, double *dx, double eps,double stoptol, double functol, int retval) {
	int rcode;
	double nrm,nrmnx,relfit;
	rcode = 0;


	if (retval == 3) {
		rcode = 4;
		return rcode;
	}
	if (retval == 15) {
		rcode = 15;
		return rcode;
	}

	nrm = l2norm(jac, N);
	nrmnx = nrm / (double) N;

	if (fabs(fo) < eps) {
		relfit = fabs(fx - fo);
	}
	else {
		relfit = fabs((fx - fo)/fo);
	}

	if (nrmnx < stoptol) {
		return 1; // Successful Convergence
	} else if (relfit < functol) {
		return 6; // Relative fit less than function tolerance
	}

	return rcode;
}

int stopcheck_mt(double fx, int N, double *xc, double *xf, double *jac, double *dx, double fsval, double gtol, double stol, int retval) {
	int rcode, i;
	double num, den;
	double stop0;
	double *scheck;

	rcode = 0;

	if (retval == 3) {
		rcode = 4;
		return rcode;
	}
	if (retval == 15) {
		rcode = 15;
		return rcode;
	}
	
	scheck = (double*)malloc(sizeof(double)*N);

	if (fabs(fx) > fabs(fsval)) {
		den = fabs(fx);
	}
	else {
		den = fabs(fsval);
	}
	for (i = 0; i < N; ++i) {
		if (fabs(xf[i]) > 1.0 / fabs(dx[i])) {
			num = fabs(xf[i]);
		}
		else {
			num = 1.0 / fabs(dx[i]);
		}
		scheck[i] = fabs(jac[i]) * num / den;
	}

	stop0 = array_max_abs(scheck, N);

	if (stop0 <= gtol) {
		rcode = 1;
	}
	else {
		for (i = 0; i < N; ++i) {
			if (fabs(xf[i]) > 1.0 / fabs(dx[i])) {
				den = fabs(xf[i]);
			}
			else {
				den = 1.0 / fabs(dx[i]);
			}
			num = fabs(xf[i] - xc[i]);
			scheck[i] = num / den;
		}
		stop0 = array_max_abs(scheck, N);
		if (stop0 <= stol) {
			rcode = 2;
		}
	}

	free(scheck);
	return rcode;
}

int cstep(double *stx,double *fx,double *dx,double *sty,double *fy,double *dy,double *stp,double *fp,double *dp,int *brackt,
                       double  stpmin,double stpmax) {
	int info,bound;
	double gamma, p, p66, q, r, s, sgnd, stpc, stpf, stpq, theta;

	info = 0;
	p66 = 0.66;

	if ((*brackt == 1 && (*stp <= pmin(*stx, *sty) || *stp >= pmax(*stx, *sty))) || *dx*(*stp - *stx) >= 0.0 || stpmax < stpmin) {
		return info;
	}

	sgnd = (*dp) * (*dx / fabs(*dx));

	/*
	 First case. A higher function value. (fp > fx)
     The minimum is bracketed. If the cubic step is closer
	 to stx than the quadratic step, the cubic step is taken,
     else the average of the cubic and quadratic steps is taken.
	*/

	if (*fp > *fx) {
		info = 1;
		bound = 1;
		theta = 3 * (*fx - *fp) / (*stp - *stx) + *dx + *dp;
		s = pmax(fabs(theta), fabs(*dx));
		s = pmax(s, fabs(*dp));
		gamma = s*sqrt((theta / s)*(theta/s) - (*dx / s)*(*dp / s));
		if (*stp < *stx) {
			gamma = -gamma;
		}
		p = (gamma - *dx) + theta;
		q = ((gamma - *dx) + gamma) + *dp;
		r = p / q;
		stpc = *stx + r*(*stp - *stx);
		stpq = *stx + ((*dx / ((*fx - *fp) / (*stp - *stx) + *dx)) / 2)*(*stp - *stx);

		if (fabs(stpc - *stx) < fabs(stpq - *stx)) {
			stpf = stpc;
		} else {
			stpf = stpc + (stpq - stpc) / 2;
		}
		*brackt = 1;
	} else if (sgnd < 0.0) {
		/*
		Second case. A lower function value and derivatives of
	    opposite sign. The minimum is bracketed. If the cubic
	    step is closer to stx than the quadratic (secant) step, 
	    the cubic step is taken, else the quadratic step is taken.
		*/

		info = 2;
		bound = 0;
		theta = 3 * (*fx - *fp) / (*stp - *stx) + *dx + *dp;
		s = pmax(fabs(theta), fabs(*dx));
		s = pmax(s, fabs(*dp));
		gamma = s*sqrt((theta / s)*(theta / s) - (*dx / s)*(*dp / s));
		if (*stp > *stx) {
			gamma = -gamma;
		}
		p = (gamma - *dp) + theta;
		q = ((gamma - *dp) + gamma) + *dx;
		r = p / q;
		stpc = *stp + r*(*stx - *stp);
		stpq = *stp + (*dp / (*dp - *dx))*(*stx - *stp);
		if (fabs(stpc - *stp) > fabs(stpq - *stp)) {
			stpf = stpc;
		} else {
			stpf = stpq;
		}
		*brackt = 1;
	}
	else if (fabs(*dp) < fabs(*dx)) {
		/*
     Third case. A lower function value, derivatives of the
     same sign, and the magnitude of the derivative decreases.
     The cubic step is only used if the cubic tends to infinity 
     in the direction of the step or if the minimum of the cubic
     is beyond stp. Otherwise the cubic step is defined to be 
     either stpmin or stpmax. The quadratic (secant) step is also 
     computed and if the minimum is bracketed then the the step 
     closest to stx is taken, else the step farthest away is taken.
		*/

		info = 3;
		bound = 1;
		theta = 3 * (*fx - *fp) / (*stp - *stx) + *dx + *dp;
		s = pmax(fabs(theta), fabs(*dx));
		s = pmax(s, fabs(*dp));

		/*
		 The case gamma = 0 only arises if the cubic does not tend
         to infinity in the direction of the step.
		*/

		gamma = s*sqrt(pmax(0., (theta / s)*(theta / s) - (*dx / s)*(*dp / s)));
		if (*stp > *stx) {
			gamma = -gamma;
		}
		p = (gamma - *dp) + theta;
		q = (gamma + (*dx - *dp)) + gamma;
		r = p / q;

		if (r < 0.0 && gamma != 0.0) {
			stpc = *stp + r*(*stx - *stp);
		} else if (*stp > *stx) {
			stpc = stpmax;
		} else {
			stpc = stpmin;
		}
		stpq = *stp + (*dp / (*dp - *dx))*(*stx - *stp);
		if (*brackt == 1) {
			if (fabs(*stp - stpc) < fabs(*stp - stpq)) {
				stpf = stpc;
			} else {
				stpf = stpq;
			}
		} else {
			if (fabs(*stp - stpc) > fabs(*stp - stpq)) {
				stpf = stpc;
			} else {
				stpf = stpq;
			}
		}
	}
	else {
		/*
     Fourth case. A lower function value, derivatives of the
     same sign, and the magnitude of the derivative does
     not decrease. If the minimum is not bracketed, the step
     is either stpmin or stpmax, else the cubic step is taken.

		*/
		info = 4;
		bound = 0;
		if (*brackt == 1) {
			theta = 3 * (*fp - *fy) / (*sty - *stp) + *dy + *dp;
			s = pmax(fabs(theta), fabs(*dy));
			s = pmax(s,fabs(*dp));
			gamma = s*sqrt((theta / s)*(theta / s) - (*dy / s)*(*dp / s));
			if (*stp > *sty) {
				gamma = -gamma;
			}
			p = (gamma - *dp) + theta;
			q = ((gamma - *dp) + gamma) + *dy;
			r = p / q;
			stpc = *stp + r*(*sty - *stp);
			stpf = stpc;
		}
		else if (*stp > *stx) {
			stpf = stpmax;
		}
		else {
			stpf = stpmin;
		}
		
	}
	/*
	Update the interval of uncertainty. This update does not
    depend on the new step or the case analysis above.
	*/

	if (*fp > *fx) {
		*sty = *stp;
		*fy = *fp;
		*dy = *dp;
	} else { 
		if (sgnd < 0.0) {
			*sty = *stx;
			*fy = *fx;
			*dy = *dx;
		}
		*stx = *stp;
		*fx = *fp;
		*dx = *dp;
	}
	/*
	Compute the new step and safeguard it.
	*/

	stpf = pmin(stpmax, stpf);
	stpf = pmax(stpmin, stpf);
	*stp = stpf;
	if (*brackt == 1 && bound == 1) {
		if (*sty > *stx) {
			*stp = pmin(*stx + p66*(*sty - *stx), *stp);
		} else {
			*stp = pmax(*stx + p66*(*sty - *stx), *stp);
		}
	}


	return info;
}

int cvsrch(custom_function *funcpt, custom_gradient *funcgrad, double *x, double *f, double *g, double *stp, double *s, int N, double *dx, double maxstep,
	int MAXITER,double eps2,double ftol, double gtol, double xtol) {
	int info,i,siter,nfev;
	int infoc, j, brackt, stage1;
	double dg, dgm, dginit, dgtest, dgx, dgxm, dgy, dgym, finit, ftest1, fm, fx, fxm, fy, fym, p5, p66, stx, sty,
		stmin, stmax, width, width1, xtrapf;
	double nlen,den,rell,stepmin;
	double *rcheck,*wa;

	rcheck = (double*)malloc(sizeof(double)*N);
	wa = (double*)malloc(sizeof(double)*N);

	nlen = 0.0;
	p5 = 0.5;
	p66 = 0.66;
	xtrapf = 4.0;
	info = 0;
	infoc = 1;
	siter = MAXITER;


	if (N <= 0 || *stp <= 0.0 || ftol < 0.0 || gtol < 0.0 || xtol < 0.0) {
		return info;
	}


	for (i = 0; i < N; ++i) {
		nlen += dx[i] * s[i] * dx[i] * s[i];
	}
	nlen = sqrt(nlen);

	for (i = 0; i < N; ++i) {
		if (fabs(x[i]) > 1.0 / fabs(dx[i])) {
			den = fabs(x[i]);
		}
		else {
			den = 1.0 / fabs(dx[i]);
		}
		rcheck[i] = s[i] / den;
	}

	rell = array_max_abs(rcheck, N);

	stepmin = ftol / rell;

	dginit = 0.0;

	for (j = 0; j < N; ++j) {
		dginit += g[j] * s[j];
	}

	if (dginit >= 0.0) {
		return info;
	}

	brackt = 0;
	stage1 = 1;
	finit = *f;
	nfev = 0;
	dgtest = ftol*dginit;
	width = maxstep - stepmin;
	width1 = width / 0.5;

	for (j = 0; j < N; ++j) {
		wa[j] = x[j];
	}

	/*     The variables stx, fx, dgx contain the values of the step, 
     function, and directional derivative at the best step.
     The variables sty, fy, dgy contain the value of the step,
     function, and derivative at the other endpoint of
     the interval of uncertainty.
     The variables stp, f, dg contain the values of the step,
     function, and derivative at the current step.
	*/

	stx = 0.0;
	fx = finit;
	dgx = dginit;
	sty = 0.0;
	fy = finit;
	dgy = dginit;

	// Iteration

	while (info == 0) {

		if (brackt == 1) {
			stmin = pmin(stx, sty);
			stmax = pmax(stx, sty);
		} else {
			stmin = stx;
			stmax = *stp + xtrapf*(*stp - stx);
		}

		*stp = pmax(*stp, stepmin);
		*stp = pmin(*stp, maxstep);

		if ((brackt == 1 && (*stp <= stmin || *stp >= stmax)) || nfev >= siter - 1 || infoc == 0 || (brackt == 1 && (stmax - stmin) <= xtol*stmax)) {
			*stp = stx;
		}

		for (j = 0; j < N; ++j) {
			x[j] = wa[j] + *stp * s[j];
		}

		*f = FUNCPT_EVAL(funcpt,x, N);
		if (*f >= DBL_MAX || *f <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			free(rcheck);
			free(wa);
			return 15;
		}
		if (*f != *f) {
			printf("Program Exiting as the function returns NaN");
			free(rcheck);
			free(wa);
			return 15;
		}
		grad_cd(funcpt,funcgrad, x, N, dx, eps2,g);
		nfev++;


		dg = 0.0;
		for (j = 0; j < N; ++j) {
			dg = dg + g[j]*s[j];
		}	
		ftest1 = finit + *stp * dgtest;

		//       Test for convergence.

		if ((brackt == 1 && (*stp <= stmin || *stp >= stmax)) || infoc == 0) {
			info = 6;
		}

		if (*stp == maxstep && *f <= ftest1 && dg <= dgtest) {
			info = 5;
		}

		if (*stp == stepmin && (*f > ftest1 || dg >= dgtest)) {
			info = 4;
		}

		if (nfev >= siter) {
			info = 3;
		}

		if (brackt == 1 && ((stmax - stmin) <= xtol*stmax)) {
			info = 2;
		}

		if (*f <= ftest1 && fabs(dg) <= gtol*(-dginit)) {
			info = 1;
		}

		if (stage1 == 1 && *f <= ftest1 && dg >= pmin(ftol, gtol)*dginit) {
			stage1 = 0;
		}


		/*
		A modified function is used to predict the step only if
        we have not obtained a step for which the modified
        function has a nonpositive function value and nonnegative
        derivative, and if a lower function value has been
        obtained but the decrease is not sufficient.
		*/
		if (stage1 == 1 && *f <= fx && *f > ftest1) {

			fm = *f - *stp*dgtest;
			fxm = fx - stx*dgtest;
			fym = fy - sty*dgtest;
			dgm = dg - dgtest;
			dgxm = dgx - dgtest;
			dgym = dgy - dgtest;

			infoc = cstep(&stx, &fxm, &dgxm, &sty, &fym, &dgym, stp, &fm, &dgm, &brackt, stmin, stmax);

			fx = fxm + stx*dgtest;
			fy = fym + sty*dgtest;
			dgx = dgxm + dgtest;
			dgy = dgym + dgtest;
		} else {
			
			infoc = cstep(&stx, &fx, &dgx, &sty, &fy, &dgy, stp, f, &dg, &brackt, stmin, stmax);
		}

		if (brackt == 1) {
			if (fabs(sty - stx) >= p66*width1) {
				*stp = stx + p5*(sty - stx);
			}
			width1 = width;
			width = fabs(sty - stx);
		}

	}

	free(rcheck);
	free(wa);

	return info;
}

int lnsrchmt(custom_function *funcpt, custom_gradient *funcgrad, double *xi, double *f, double *jac, double *alpha, double *p, int N, double *dx, double maxstep, int MAXITER,
		double eps2,double ftol, double gtol, double xtol, double *x) {
	int i,retval,info;

	/*
	* Return Codes
	*
	* Codes 1,2,3 denote possible success.
	* Codes 0 and 4 denote failure.
	*
	* 1 - df(x)/dx <= gtol achieved so xf may be the local minima.
	* 2 - Distance between the last two steps is less than stol or |xf - xi| <= stol so xf may be the local minima.
	* 3 - Global Step failed to locate a lower point than xi so xi may be the local minima.
	* 4 - Iteration Limit exceeded. Convergence not achieved.
	*
	*/

	//f = funcpt(xi, N);

	// Important - All argument values are modified in this algorithm 

	for (i = 0; i < N; ++i) {
		x[i] = xi[i];
	}

	info = cvsrch(funcpt,funcgrad, x,f, jac, alpha, p, N, dx, maxstep,MAXITER,eps2,
		ftol, gtol, xtol);
	
	retval = info;
	return retval;

}
