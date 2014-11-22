/*
 * optimize.c
 *
 *  Copyright (c) 2014, Rafat Hussain
 *	License : BSD 3-Clause
 *	See COPYRIGHT for more details
 */
 


#include "optimc.h"

opt_object opt_init(int N) {
	opt_object obj = NULL;
	int i;
	double meps;

	obj = (opt_object) malloc (sizeof(struct opt_set) + sizeof(double)* (N-1));

	meps = macheps();
	obj->eps = meps;
	obj->xtol = meps;
	obj->gtol = pow(meps,(double)1.0/3.0);
	obj->ftol = obj->gtol * obj->gtol;
	obj->stol = obj->ftol;
	obj->N = N;
	obj->MaxIter = 200*N;
	obj->retval = 0;

	if (obj->MaxIter < 1000) {
		obj->MaxIter = 1000;
	}

	obj->Iter = 0;
	obj->Method = 0;
	obj->maxstep = -1.0;
	strcpy(obj->MethodName,"Nelder-Mead");
	obj->objfunc = 0.0;
	for (i = 0; i < N;++i) {
		obj->xopt[i] = 0.0;
	}

	return obj;
}

void optsummary(opt_object obj) {
	int i;
	printf("\n Return Value : %d \n",obj->retval);
	printf("Method : %d %s \n",obj->Method,obj->MethodName);
	printf("Iterations : %d \n",obj->Iter);
	printf("Function Minimized At : [ ");
	for(i = 0; i < obj->N;++i) {
		printf("%g ",obj->xopt[i]);
	}
	printf(" ] \n");
	printf("Function Value : %g \n \n",obj->objfunc);

}

void setMaxIter(opt_object obj,int MaxIter) {
	obj->MaxIter = MaxIter;
}

void setMaxStep(opt_object obj, double maxstep) {
	obj->maxstep = maxstep;
}

void setTOL(opt_object obj,double gtol,double stol,double ftol,double xtol) {
	obj->gtol = gtol;
	obj->stol = stol;
	obj->ftol = ftol;
	obj->xtol = xtol;
}

nls_object nls_init(int M,int N) {
	nls_object obj = NULL;
	int i;
	double meps;

	obj = (nls_object) malloc (sizeof(struct nls_set) + sizeof(double)* (N-1));

	meps = macheps();
	obj->eps = meps;
	obj->epsfcn = meps;
	obj->factor = 100.; //Takes values between 0.1 and 100.0
	obj->xtol = meps;
	obj->gtol = pow(meps,(double)1.0/3.0);
	obj->ftol = obj->gtol * obj->gtol;
	obj->M = M;
	obj->N = N;
	obj->MaxIter = 1000;
	obj->Maxfev = 1000;
	obj->nfev = 0;
	obj->njev = 0;
	obj->ldfjac = M;
	obj->mode = 1;

	if (obj->MaxIter < 1000) {
		obj->MaxIter = 1000;
	}

	obj->Iter = 0;
	for (i = 0; i < N;++i) {
		obj->xopt[i] = 0.0;
	}

	return obj;
}

void setnlsTOL(nls_object obj,double gtol,double ftol,double xtol) {
	obj->gtol = gtol;
	obj->ftol = ftol;
	obj->xtol = xtol;
}

int fminsearch(custom_function *funcpt,int N,double *xi,double *xf) {
	int i,retval,MAXITER,niter;
	double fsval,eps;
	double *dx;

	dx = (double*) malloc(sizeof(double) * N);

	fsval = 1.0;
	MAXITER = 200*N;
	niter = 0;
	eps = macheps(); // Use macheps program

	for(i = 0; i < N;++i) {
		dx[i] = 1.0;
	}

	retval = nel_min(funcpt,xi,N,dx,fsval,MAXITER,&niter,eps,xf);

	//printf("Iterations %d \n", niter);

	free(dx);
	return retval;
}

static int mvalue(int N) {
	int mval;

	if (N <= 10) {
		mval = N;
	} else if (N > 10 && N <= 20) {
		mval = 10;
	} else if (N > 20 && N <= 200) {
		mval = 15;
	} else if ( N > 200) {
		mval = 20;
	}

	return mval;
}

int fminunc(custom_function *funcpt, custom_gradient *funcgrad, int N, double *xi,double maxstep, int method,double *xf) {
	int i,retval,MAXITER,niter,m;
	double fsval,eps,gtol,stol,ftol,xtol,delta;
	double *dx;

	dx = (double*) malloc(sizeof(double) * N);
		/*
		 * Method 0 - Nelder-Mead
		 * Method 1 - Newton Line Search
		 * Method 2 - Newton Trust Region - Hook Step
		 * Method 3 - Newton Trust Region - Double Dog-Leg
		 * Method 4 - Conjugate Gradient
		 * Method 5 - BFGS
		 * Method 6 - Limited Memory BFGS
		 * Method 7 - BFGS More-Thuente Line Search
		 */
	/*
	 * Return Codes
	 *
	 * Code 1 denotes probable success.
	 * Codes 2,3 and 6 denote possible success.
	 * Codes 0, 4 and 15 denote failure. [Code 4 my also occur in cases where the minima is realized but
	 * convergence is not achieved due to tolerance values being too small. It is recommended that you use
	 * another method if you suspect that such a scenario has occured.]
	 *
	 * 0 - Input Error
	 * 1 - df(x)/dx <= gtol achieved so current point may be the local minima.
	 * 2 - Distance between the last two steps is less than stol or |xf - xi| <= stol so the point may be the local minima.
	 * 3 - Global Step failed to locate a lower point than the current point so it may be the local minima.
	 * 4 - Iteration Limit exceeded. Convergence probably not achieved.
	 * 6 - Function value drops below ftol (relative functional tolerance).
	 * 15 - Overflow occurs. Try a different method.
	 *
	 */
	fsval = 1.0;
	MAXITER = 200*N;
	niter = 0;
	delta = -1.0; // Trust Region default
	eps = macheps(); // Use macheps program

	for(i = 0; i < N;++i) {
		dx[i] = 1.0;
	}

	if (method == 0) {
		retval = nel_min(funcpt,xi,N,dx,fsval,MAXITER,&niter,eps,xf);
	} else if (method == 1) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_func(funcpt,funcgrad,xi,N,dx,fsval,maxstep,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 2) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_trust(funcpt,funcgrad,xi,N,dx,fsval,delta,0,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 3) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_trust(funcpt,funcgrad,xi,N,dx,fsval,delta,1,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 4) {
		gtol = pow(eps,1.0/3.0);
		ftol = gtol * gtol;
		xtol = eps;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = conjgrad_min_lin(funcpt,funcgrad,xi,N,dx,maxstep,MAXITER,&niter,eps,gtol,ftol,xtol,xf);
	} else if (method == 5) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = bfgs_min(funcpt,funcgrad,xi,N,dx,fsval,maxstep,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 6) {
		gtol = pow(eps,1.0/3.0);
		ftol = gtol * gtol;
		xtol = eps;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		m = mvalue(N);
		retval = bfgs_l_min(funcpt,funcgrad,xi,N,m,dx,fsval,maxstep,MAXITER,&niter,eps,gtol,ftol,xtol,xf);
	}
	else if (method == 7) {
		gtol = pow(eps, 1.0 / 3.0);
		ftol = gtol * gtol;
		xtol = eps;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		m = mvalue(N);
		retval = bfgs_min2(funcpt, funcgrad, xi, N, m, dx, fsval, maxstep, MAXITER, &niter, eps, gtol, ftol, xtol, xf);
	}
	else {
		printf("Method Value should be one of 0,1,2,3,4,5,6 or 7. See Documentation. \n");
		exit(1);
	}



	free(dx);
	return retval;
}

double fminbnd(custom_funcuni *funcuni,double a, double b) {
	double x,t,eps;

	t = 1e-012;
	eps = macheps();

	brent_local_min(funcuni,a,b,t,eps,&x);

	return x;
}

int fminnewt(custom_function *funcpt, custom_gradient *funcgrad, int N, double *xi,
	     double delta,double *dx,double fsval,double maxstep, int method,double *xf) {
	int retval;
	int MAXITER,niter;
	double eps,gtol,stol;
	/*

	 * Method 1 - Newton Line Search
	 * Method 2 - Newton Trust Region - Hook Step
	 * Method 3 - Newton Trust Region - Double Dog-Leg
	 *
	 * Default Values :
	 *
	 * fsval = 1.0
	 * delta = -1.0
	 * dx = {1.0,1.0,...} - 1XN vector

	 */
	/*
	 * Return Codes
	 *
	 * Code 1 denotes probable success.
	 * Codes 2,3 and 6 denote possible success.
	 * Codes 0, 4 and 15 denote failure. [Code 4 my also occur in cases where the minima is realized but
	 * convergence is not achieved due to tolerance values being too small. It is recommended that you use
	 * another method if you suspect that such a scenario has occured.]
	 *
	 * 0 - Input Error
	 * 1 - df(x)/dx <= gtol achieved so current point may be the local minima.
	 * 2 - Distance between the last two steps is less than stol or |xf - xi| <= stol so the point may be the local minima.
	 * 3 - Global Step failed to locate a lower point than the current point so it may be the local minima.
	 * 4 - Iteration Limit exceeded. Convergence probably not achieved.
	 * 6 - Function value drops below ftol (relative functional tolerance).
	 * 15 - Overflow occurs. Try a different method.
	 *
	 */

	MAXITER = 200*N;
	niter = 0;

	eps = macheps(); // Use macheps program

	if (method == 1) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_func(funcpt,funcgrad,xi,N,dx,fsval,maxstep,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 2) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_trust(funcpt,funcgrad,xi,N,dx,fsval,delta,0,MAXITER,&niter,eps,gtol,stol,xf);
	} else if (method == 3) {
		gtol = pow(eps,1.0/3.0);
		stol = gtol * gtol;
		if (MAXITER < 1000) {
			MAXITER = 1000;
		}
		retval = newton_min_trust(funcpt,funcgrad,xi,N,dx,fsval,delta,1,MAXITER,&niter,eps,gtol,stol,xf);
	} else {
		printf("Method Value should be one of 1,2 or 3. See Documentation. \n");
		exit(1);
	}


	return retval;
}

void optimize(opt_object obj, custom_function *funcpt, custom_gradient *funcgrad, int N, double *xi,
		int method) {
	int i,m;
	double fsval,delta;
	double *dx;

	dx = (double*) malloc(sizeof(double) * N);
		/*
		 * Method 0 - Nelder-Mead
		 * Method 1 - Newton Line Search
		 * Method 2 - Newton Trust Region - Hook Step
		 * Method 3 - Newton Trust Region - Double Dog-Leg
		 * Method 4 - Conjugate Gradient
		 * Method 5 - BFGS
		 * Method 6 - Limited Memory BFGS
		 * Method 7 - BFGS More-Thuente Line Search
		 */
	/*
	 * Return Codes
	 *
	 * Code 1 denotes probable success.
	 * Codes 2,3 and 6 denote possible success.
	 * Codes 0, 4 and 15 denote failure. [Code 4 my also occur in cases where the minima is realized but
	 * convergence is not achieved due to tolerance values being too small. It is recommended that you use
	 * another method if you suspect that such a scenario has occurred.]
	 *
	 * 0 - Input Error
	 * 1 - df(x)/dx <= gtol achieved so current point may be the local minima.
	 * 2 - Distance between the last two steps is less than stol or |xf - xi| <= stol so the point may be the local minima.
	 * 3 - Global Step failed to locate a lower point than the current point so it may be the local minima.
	 * 4 - Iteration Limit exceeded. Convergence probably not achieved.
	 * 6 - Function value drops below ftol (relative functional tolerance).
	 * 15 - Overflow occurs. Try a different method.
	 *
	 */
	fsval = 1.0;
	delta = -1.0; // Trust Region default


	for(i = 0; i < N;++i) {
		dx[i] = 1.0;
	}
	obj->Method = method;
	obj->Iter = 0;

	if (obj->N != N) {
		printf("The Object is initialized for a problem of size %d \n",obj->N);
		printf("Please Reinitialize the object for the new size %d using opt_init command \n",N);
		exit(1);
	}

	if (method == 0) {
		strcpy(obj->MethodName,"Nelder-Mead");
		obj->retval = nel_min(funcpt,xi,obj->N,dx,fsval,obj->MaxIter,&obj->Iter,obj->eps,obj->xopt);
	} else if (method == 1) {
		strcpy(obj->MethodName,"Newton Line Search");
		obj->retval = newton_min_func(funcpt,funcgrad,xi,obj->N,dx,fsval,obj->maxstep,obj->MaxIter,&obj->Iter,obj->eps,obj->gtol,
				obj->stol,obj->xopt);
	} else if (method == 2) {
		strcpy(obj->MethodName,"Newton Trust Region - Hook Step");

		obj->retval = newton_min_trust(funcpt,funcgrad,xi,obj->N,dx,fsval,delta,0,obj->MaxIter,&obj->Iter,obj->eps,
				obj->gtol,obj->stol,obj->xopt);
	} else if (method == 3) {
		strcpy(obj->MethodName,"Newton Trust Region - Double Dog-Leg");

		obj->retval = newton_min_trust(funcpt,funcgrad,xi,obj->N,dx,fsval,delta,1,obj->MaxIter,&obj->Iter,obj->eps,
				obj->gtol,obj->stol,obj->xopt);
	} else if (method == 4) {
		strcpy(obj->MethodName,"Conjugate Gradient");

		obj->retval = conjgrad_min_lin(funcpt, funcgrad, xi, obj->N, dx, obj->maxstep, obj->MaxIter, &obj->Iter, obj->eps, obj->gtol,
				obj->ftol,obj->xtol,obj->xopt);
	} else if (method == 5) {
		strcpy(obj->MethodName,"BFGS");

		obj->retval = bfgs_min(funcpt, funcgrad, xi, obj->N, dx, fsval, obj->maxstep, obj->MaxIter, &obj->Iter, obj->eps, obj->gtol,
				obj->stol,obj->xopt);
	} else if (method == 6) {
		strcpy(obj->MethodName,"Limited Memory BFGS");
		m = mvalue(N);
		obj->retval = bfgs_l_min(funcpt, funcgrad, xi, obj->N, m, dx, fsval, obj->maxstep, obj->MaxIter, &obj->Iter, obj->eps,
				obj->gtol,obj->ftol,obj->xtol,obj->xopt);
	}
	else if (method == 7) {
		strcpy(obj->MethodName, "BFGS More-Thuente Line Search");
		m = mvalue(N);
		obj->retval = bfgs_min2(funcpt, funcgrad, xi, obj->N, m, dx, fsval, obj->maxstep, obj->MaxIter, &obj->Iter, obj->eps,
			obj->gtol, obj->ftol, obj->xtol, obj->xopt);
	}
	else {
		strcpy(obj->MethodName,"NULL");
		printf("Method Value should be one of 0,1,2,3,4,5 or 6. See Documentation. \n");
		exit(1);
	}

	obj->objfunc = FUNCPT_EVAL(funcpt,obj->xopt,obj->N);

	free(dx);
}

void free_opt(opt_object object) {

	free(object);

}

int levmar(custom_funcmult *funcmult, custom_jacobian *jacobian,
		double *xi,int M, int N,double *xf) {
	int info,i;
	double *fvec,*fjac,*diag,*qtf;
	int ldfjac,mode,nfev,njev,nprint,maxfev;
	int *ipvt;
	double factor,eps,epsfcn,ftol,gtol,xtol;

	fvec = (double*) malloc(sizeof(double) *M);
	fjac = (double*) malloc(sizeof(double) *M * N);
	diag = (double*) malloc(sizeof(double) *N);
	qtf = (double*) malloc(sizeof(double) *N);
	ipvt = (int*) malloc(sizeof(int) *N);

	eps = macheps();
	epsfcn = eps;
	gtol = pow(eps,1.0/3.0);
	ftol = gtol * gtol;
	xtol = eps;
	factor = 100.;
	ldfjac = M;
	nfev = njev = nprint = 0;
	mode = 1;
	maxfev = 1000;

	for(i = 0; i < N;++i) {
		xf[i] = xi[i];
	}

	if (jacobian == NULL) {
		info = lmdif(funcmult,xf,M,N,fvec,fjac,ldfjac,maxfev,diag,mode,factor,nprint,eps,epsfcn,ftol,gtol,
				xtol,&nfev,&njev,ipvt,qtf);
	} else {
		info = lmder(funcmult,jacobian,xf,M,N,fvec,fjac,ldfjac,maxfev,diag,mode,factor,nprint,eps,ftol,gtol,
				xtol,&nfev,&njev,ipvt,qtf);
	}


	free(fvec);
	free(fjac);
	free(diag);
	free(qtf);
	free(ipvt);

	return info;
}

void nls(nls_object obj, custom_funcmult *funcmult, custom_jacobian *jacobian,
		double *xi) {

	int info,i,M,N,nprint;
	double *fvec,*fjac,*diag,*qtf;
	int *ipvt;

	M = obj->M;
	N = obj-> N;
	nprint = 0;

	fvec = (double*) malloc(sizeof(double) *M);
	fjac = (double*) malloc(sizeof(double) *M * N);
	diag = (double*) malloc(sizeof(double) *N);
	qtf = (double*) malloc(sizeof(double) *N);
	ipvt = (int*) malloc(sizeof(int) *N);

	for(i = 0; i < N;++i) {
		obj->xopt[i] = xi[i];
	}

	if (jacobian == NULL) {
		info = lmdif(funcmult,obj->xopt,obj->M,obj->N,fvec,fjac,obj->ldfjac,obj->Maxfev,diag,obj->mode,
				obj->factor,nprint,obj->eps,obj->epsfcn,obj->ftol,obj->gtol,obj->xtol,&obj->nfev,&obj->njev,ipvt,qtf);
	} else {
		info = lmder(funcmult,jacobian,obj->xopt,obj->M,obj->N,fvec,fjac,obj->ldfjac,obj->Maxfev,diag,obj->mode,
				obj->factor,nprint,obj->eps,obj->ftol,obj->gtol,obj->xtol,&obj->nfev,&obj->njev,ipvt,qtf);
	}
	obj->retval = info;
	obj->Iter = obj->nfev;

	free(fvec);
	free(fjac);
	free(diag);
	free(qtf);
	free(ipvt);

}

void nls_scale(nls_object obj, custom_funcmult *funcmult, custom_jacobian *jacobian,
		double *diag,double *xi) {

	int info,i,M,N,nprint;
	double *fvec,*fjac,*qtf;
	int *ipvt;

	M = obj->M;
	N = obj-> N;
	obj->mode = 2;
	nprint = 0;

	fvec = (double*) malloc(sizeof(double) *M);
	fjac = (double*) malloc(sizeof(double) *M * N);
	qtf = (double*) malloc(sizeof(double) *N);
	ipvt = (int*) malloc(sizeof(int) *N);

	for(i = 0; i < N;++i) {
		obj->xopt[i] = xi[i];
	}

	if (jacobian == NULL) {
		info = lmdif(funcmult,obj->xopt,obj->M,obj->N,fvec,fjac,obj->ldfjac,obj->Maxfev,diag,obj->mode,
				obj->factor,nprint,obj->eps,obj->epsfcn,obj->ftol,obj->gtol,obj->xtol,&obj->nfev,&obj->njev,ipvt,qtf);
	} else {
		info = lmder(funcmult,jacobian,obj->xopt,obj->M,obj->N,fvec,fjac,obj->ldfjac,obj->Maxfev,diag,obj->mode,
				obj->factor,nprint,obj->eps,obj->ftol,obj->gtol,obj->xtol,&obj->nfev,&obj->njev,ipvt,qtf);
	}
	obj->retval = info;
	obj->Iter = obj->nfev;

	free(fvec);
	free(fjac);
	free(qtf);
	free(ipvt);

}

void nlssummary(nls_object obj) {
	int i;
	printf("\n\n");
	printf("Return Value : %d \n",obj->retval);
	if (obj->retval == 0) {
		printf("Termination Message : improper input parameters. \n");
	}
	if (obj->retval == 1) {
		printf("Termination Message : both actual and predicted relative reductions "
				"in the sum of squares are at most ftol. \n");
	}
	if (obj->retval == 2) {
		printf("Termination Message : relative error between two consecutive iterates "
				"is at most xtol. \n");
	}
	if (obj->retval == 3) {
		printf("Termination Message : conditions for info = 1 and info = 2 both hold. \n");
	}
	if (obj->retval == 4) {
		printf("Termination Message : the cosine of the angle between fvec and any "
				"column of the jacobian is at most gtol in absolute value. \n");
	}
	if (obj->retval == 5) {
		printf("Termination Message : number of calls to fcn has reached or "
				"exceeded maxfev. \n");
	}
	if (obj->retval == 6) {
		printf("Termination Message : ftol is too small. no further reduction in "
				"the sum of squares is possible. \n");
	}
	if (obj->retval == 7) {
		printf("Termination Message : xtol is too small. no further improvement in"
				"the approximate solution x is possible. \n");
	}
	if (obj->retval == 8) {
		printf("Termination Message : gtol is too small. fvec is orthogonal to the"
				"columns of the jacobian to machine precision. \n");
	}
	printf("Iterations : %d \n",obj->Iter);
	printf("Function Minimized At : [ ");
	for(i = 0; i < obj->N;++i) {
		printf("%g ",obj->xopt[i]);
	}
	printf(" ] \n");
	//printf("Function Value : %g \n \n",obj->objfunc);

}

void free_nls(nls_object object) {

	free(object);

}

