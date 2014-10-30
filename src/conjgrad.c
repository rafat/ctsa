#include "conjgrad.h"

static int richolu(double *A,int N, int stride, double *U22) {
	int sc;
	int j,i,u,w;
	double u11;
	
	if (N == 1) {
		if (A[0] > 0) {
			A[0] = sqrt(A[0]);
			return 0;
		} else {
			return -1;
		}
	} else {
		if (A[0] < 0) {
			return -1;
		}
		u11 = sqrt(A[0]);
		A[0] = u11;
		for (j = 1; j < N;++j) {
			if (A[j] != 0) {
				A[j] /= u11;
			}
		}
		mmult(A+1,A+1,U22,N-1,1,N-1);
		for (i = 0; i < N-1; ++i) {
			u = stride + 1 + i * stride;
			w = i * (N-1);
			for(j = i; j < N-1;j++) {
				if (A[j + u] != 0) {
					A[j + u] -= U22[j + w];
				}
			}
		}
		
		sc = richolu(A+stride+1,N-1,stride,U22);
		if (sc == -1) {
			return -1;
		}
		
	}
	
	return sc;
	
}


static int icholu(double *A, int N) {
	int stride,i,j,t,sc;
	double *U22;
	U22 = (double*) malloc(sizeof(double) * N * N);
	stride = N; 
	
	sc = richolu(A,N,stride,U22);
	
	for(i=0; i < N;++i) {
		t = i *N;
		for(j=0;j < i;++j) {
			A[t+j] = 0.;
		}
	}

	free(U22);
	return sc;
	
}

int ichol(double *A, int N) {
	int sc;
	sc = icholu(A,N);
	return sc;
}

int stopcheck2(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol) {
	int rcode,i;
	double num,den;
	double stop0;
	double *scheck;
	
	rcode = 0;	
	
	scheck = (double*) malloc(sizeof(double) *N);
	
	if (fabs(fx) > fabs(fsval)) {
			den = fabs(fx);
	} else {
			den = fabs(fsval);
	}
	for(i = 0; i < N;++i) {
		if (fabs(xf[i]) > 1.0 / fabs(dx[i])) {
			num = fabs(xf[i]);
		} else {
			num = 1.0 / fabs(dx[i]);
		}
		scheck[i] = fabs(jac[i]) * num / den;
	}
	
	stop0 = array_max_abs(scheck,N);
	
	if (stop0 <= gtol) {
		rcode = 1;
	} else {
		for(i = 0; i < N;++i) {
			if (fabs(xf[i]) > 1.0 / fabs(dx[i])) {
				den = fabs(xf[i]);
			} else {
				den = 1.0 / fabs(dx[i]);
			}
			num = fabs(xf[i] - xc[i]);
			scheck[i] = num / den;
		}
		stop0 = array_max_abs(scheck,N);
		if (stop0 <= stol) {
			rcode = 2;
		}
	}
	
	free(scheck);
	return rcode;
}


int cgpr_mt(custom_function *funcpt, custom_gradient *funcgrad, double *xc, int N, double *dx,double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double ftol,double xtol,double *xf) {
	int i, rcode, retval, k, restart,gfdcode;
	int siter;
	double *xi,*temp, *rk, *pk, *jac, *jacf, *apk;
	double fsval,fxf,eps2,fo;
	double dt1, dt2,alpha;

	temp = (double*)malloc(sizeof(double)* 8);
	rk = (double*)malloc(sizeof(double)*N);
	pk = (double*)malloc(sizeof(double)*N);
	apk = (double*)malloc(sizeof(double)*N);
	jac = (double*)malloc(sizeof(double)*N);
	jacf = (double*)malloc(sizeof(double)*N);
	xi = (double*)malloc(sizeof(double)*N);

	*niter = 0;
	k = 0;
	rcode = 0;
	restart = N;
	siter = MAXITER;
	eps2 = sqrt(eps);
	gfdcode = 0;

	//xtol = 1.0e-15;
	//ftol = 1e-10;
	//gtol = 1e-05;

	// Values that may not be needed

	fsval = 1.0;
	alpha = 1.0;

	for (i = 0; i < N; ++i) {
		xi[i] = xc[i];
	}
	
	if (maxstep <= 0.0) {
		maxstep = 1000.0;
		dt1 = dt2 = 0.0;
		for (i = 0; i < N; ++i) {
			dt1 += dx[i] * dx[i];
			dt2 += dx[i] * xi[i] * dx[i] * xi[i];
		}

		dt1 = sqrt(dt1);
		dt2 = sqrt(dt2);
		if (dt1 > dt2) {
			maxstep *= dt1;
		}
		else {
			maxstep *= dt2;
		}
	}



	gfdcode = grad_fd(funcpt,funcgrad, xi, N, dx,eps2, jac);
	if (gfdcode == 15) {
		rcode = 15;
	}
	for (i = 0; i < N; ++i) {
		pk[i] = -jac[i];
		xf[i] = xi[i];
	}

	fxf = FUNCPT_EVAL(funcpt,xi, N);

	if (fxf >= DBL_MAX || fxf <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		rcode = 15;
	}
	if (fxf != fxf) {
		printf("Program Exiting as the function returns NaN");
		rcode = 15;
	}

	if (restart < N) {
		restart = N;
	}

	fo = fxf;


	while (rcode == 0 && *niter < siter) {
		*niter = *niter + 1;
		k++;

		mmult(jac, jac, temp, 1, N, 1);
		for (i = 0; i < N; ++i) {
			jacf[i] = jac[i];
		}

		retval = lnsrchmt(funcpt,funcgrad, xi,&fxf,jac,&alpha, pk, N, dx, maxstep,MAXITER,eps2, ftol, gtol, xtol, xf);
		if (retval == 100) {
			printf("The Linesearch Algorithm didn't converge");
			break;
		}
		//mdisplay(xf,1,N);
		
		//mdisplay(xf, 1, N);
		//grad_fd(funcpt, xf, N, dx, jacf);
		mmult(jac, jac, temp + 1, 1, N, 1);
		mmult(jacf, jac, temp + 3, 1, N, 1);
		temp[2] = (temp[1] - temp[3]) / temp[0]; // beta

		if (temp[2] < 0) {
			temp[2] = 0;
		}
		if (k == restart) {
			for (i = 0; i < N; ++i) {
				pk[i] = -jac[i];
			}
			k = 0;
		}
		else {
			for (i = 0; i < N; ++i) {
				pk[i] = temp[2] * pk[i] - jac[i];
			}
		}
		

		rcode = stopcheck2_mt(fxf,N,fo,jac,dx,eps,gtol,ftol,retval);
		fo = fxf;
		for (i = 0; i < N; ++i) {
			xi[i] = xf[i];
		}
	}

	if (rcode == 0 && *niter >= siter) {
		rcode = 4;
	}

	free(temp);
	free(rk);
	free(pk);
	free(jac);
	free(apk);
	free(jacf);
	free(xi);

	return rcode;
}

int conjgrad_min_lin(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double ftol,double xtol,double *xf) {
	int rcode,i;
	
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
	 * 5.- Only First Strong Wolfe Condition Achieved (Using Conjugate Gradient Method)
	 * 
	 */ 

	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	
	//rcode = cgpc(xi,N,A,b,xf);
	//rcode = cgpr(funcpt,xi,N,dx,xf);// FR
	rcode = cgpr_mt(funcpt,funcgrad,xi,N,dx,maxstep,MAXITER,niter,eps,gtol,ftol,xtol,xf);//PR+

	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	
	
	return rcode;
}
