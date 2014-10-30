#include "secant.h"

/*
 * secant.c
 *
 *  Copyright (c) 2014, Rafat Hussain
 *	License : BSD 3-Clause
 *	See COPYRIGHT for more details
 */

static void jrotate(double *A,int N,double a,double b,int i) {
	int j,r,r1;
	double c,s,den,y,w;
	r = i * N;
	r1 = r + N;
	if (a == 0.0) {
		c = 0.0;
		s = signx(b);
	} else {
		den = sqrt(a*a+b*b);
		c = a / den;
		s = b / den;
	}
	
	for (j = i; j < N;++j) {
		y = A[r+j];
		w = A[r1+j];
		A[r+j] = c*y - s*w;
		A[r1+j] = s*y + c*w;
	}
}

static void qrupdate(double *R,int N,double *u,double *v) {
	int i,j,k,l;
	
	for(i = 1; i < N;++i) {
		R[i*N+i-1] = 0;
	}
	
	k = N-1;
	while (u[k] == 0 && k > 0) {
		k--;
	}
	
	for (i = k-1; i >= 0;--i) {
		jrotate(R,N,u[i],-u[i+1],i);
		if (u[i] == 0) {
			u[i] = fabs(u[i+1]);
		} else {
			u[i] = sqrt(u[i]*u[i] + u[i+1]*u[i+1]);
		}
	}
	
	for(i = 0; i < N;++i) {
		R[i] += u[0] * v[i];
	}
	
	for(i = 0; i < k;++i) {
		j = i * N;
		l = j + N;
		jrotate(R,N,R[j+i],-R[l+i],i);
	}
}

void bfgs_naive(double *H,int N,double eps,double *xi,double *xf,double *jac,double *jacf) {
	int i,j,supd,ct;
	double sn,yn,fd,yt,jacm;
	double *sk,*yk,*temp,*temp2,*t;
	
	sk = (double*) malloc(sizeof(double) *N);
	yk = (double*) malloc(sizeof(double) *N);
	temp = (double*) malloc(sizeof(double) *1);
	temp2 = (double*) malloc(sizeof(double) *1);
	t = (double*) malloc(sizeof(double) *N);
	
	msub(xf,xi,sk,1,N);
	msub(jacf,jac,yk,1,N);
	
	sn = l2norm(sk,N);
	yn = l2norm(yk,N);
	fd = sqrt(eps);
	
	mmult(yk,sk,temp,1,N,1);
	
	if (temp[0] >= fd*sn*yn) {
		supd = 1;
		for(i = 0; i < N;++i) {
			t[i] = 0.0;
			for(j = 0; j < i;++j) {
				ct = j * N;
				t[i] += H[ct+i] * sk[j];
			}
			
			for(j = i; j < N;++j) {
				ct = i * N;
				t[i] += H[ct + j] * sk[j];
			}
			yt = fabs(yk[i] - t[i]);
			if (jac[i] > jacf[i]) {
				jacm = jac[i];
			} else {
				jacm = jacf[i];
			}
			if (yt >= fd * jacm) {
				supd = 0;
			}
		}
		if (supd == 0) {
			mmult(sk,t,temp2,1,N,1);
			for(i = 0;i < N; ++i) {
				ct = i * N;
				for(j = i; j < N;++j) {
					H[ct+j] += (yk[i]*yk[j]/temp[0]);
					H[ct+j] -= (t[i]*t[j]/temp2[0]);
				}
			}
		} 
	} 
	
	free(sk);
	free(yk);
	free(temp);
	free(temp2);
	free(t);
}

static void inithess_naive(double *H,int N,double fi,double fsval,double *dx) {
	int i,j,ct;
	double temp;
	
	if (fabs(fi) > fsval) {
		temp = fabs(fi);
	} else {
		temp = fsval;
	}
	
	for(i = 0; i < N;++i) {
		ct = i *N;
		for(j = 0; j < N;++j) {
			if (i == j) {
				H[ct+j] = temp * dx[i] * dx[i];
			} else {
				H[ct+j] = 0.0;
			}
			
		}
	}
	
}

int bfgs_min_naive(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double fsval, double maxstep, int MAXITER,
		double eps,double *xf)  {
	int rcode,iter,gfdcode;
	int i,siter,retval;
	double gtol,stol,dt1,dt2;
	double fx,num,den,stop0,fxf,eps2;
	double *jac,*hess,*scheck,*xc,*L,*step,*jacf;
	
	jac = (double*) malloc(sizeof(double) *N);
	scheck = (double*) malloc(sizeof(double) *N);
	xc = (double*) malloc(sizeof(double) *N);
	step = (double*) malloc(sizeof(double) *N);
	hess = (double*) malloc(sizeof(double) *N * N);
	L = (double*) malloc(sizeof(double) *N * N);
	jacf = (double*) malloc(sizeof(double) *N);
	
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
	
	rcode = 0;
	iter = 0;
	siter = MAXITER;
	gtol = pow(eps,1.0/3.0);
	
	stol = gtol * gtol;
	eps2 = sqrt(eps);
	gfdcode = 0;
	//set values
	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	fx = FUNCPT_EVAL(funcpt,xi,N);
	if (fx >= DBL_MAX || fx <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		rcode = 15;
	}
	if (fx != fx) {
		printf("Program Exiting as the function returns NaN");
		rcode = 15;
	}
	
	gfdcode = grad_fd(funcpt,funcgrad,xi,N,dx,eps2,jac);
	if (gfdcode == 15) {
		rcode = 15;
	}
	
	
	//maxstep = 1000.0; // Needs to be set at a much higher value proportional to l2 norm of dx

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
	
	//printf("dt1 dt2 %g \n", maxstep);

	//Check Stop0
	if (fabs(fx) > fabs(fsval)) {
			den = fabs(fx);
	} else {
			den = fabs(fsval);
	}
	for(i = 0; i < N;++i) {
		if (fabs(xi[i]) > 1.0 / fabs(dx[i])) {
			num = fabs(xi[i]);
		} else {
			num = 1.0 / fabs(dx[i]);
		}
		scheck[i] = fabs(jac[i]) * num / den;
	}
	
	stop0 = array_max_abs(scheck,N);
	
	if (stop0 <= gtol * 1e-03) {
		rcode = 1;
		for(i = 0; i < N;++i) {
			xf[i] = xi[i];
		}
		return rcode;
	}
	
	//hessian_fd(funcpt,xi,N,dx,hess);
	inithess_naive(hess,N,fx,fsval,dx);
	
	for(i = 0; i < N;++i) {
		xc[i] = xi[i];
	}
	
	while (rcode == 0 && iter < siter) {
		iter++;
		modelhess(hess,N,dx,eps,L);
		scale(jac,1,N,-1.0);
		//mdisplay(hess,N,N);
		
		linsolve_lower(L,N,jac,step);
		
		scale(jac,1,N,-1.0);
		//retval = lnsrchmod(funcpt,xc,jac,step,N,dx,maxstep,stol,xf,jacf); 
		retval = lnsrch(funcpt,xc,jac,step,N,dx,maxstep,stol,xf);
		
		//retval = swolfe(funcpt,xc,jac,step,N,dx,maxstep,stol,xf);
		
		fxf = FUNCPT_EVAL(funcpt,xf,N);
		if (fxf >= DBL_MAX || fxf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			rcode = 15;
			break;
		}
		if (fxf != fxf) {
			printf("Program Exiting as the function returns NaN");
			rcode = 15;
			break;
		}
		//printf("%d %g \n",iter,fxf);
		gfdcode = grad_fd(funcpt,funcgrad,xf,N,dx,eps2,jacf);
		if (gfdcode == 15) {
			rcode = 15;
			break;
		}
		rcode = stopcheck(fxf,N,xc,xf,jacf,dx,fsval,gtol,stol,retval);
		//hessian_fd(funcpt,xf,N,dx,hess);
		bfgs_naive(hess,N,eps,xc,xf,jac,jacf);
		for(i = 0; i < N;++i) {
			xc[i] = xf[i];
			jac[i] = jacf[i];
		}
	}
	
	if (rcode == 0 && iter >= siter) {
		rcode = 4;
	}
	/*
	for(i = 0; i < N;++i) {
		xf[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	*/
	free(jac);
	free(hess);
	free(scheck);
	free(xc);
	free(L);
	free(step);
	free(jacf);
	return rcode;
}

static void inithess_lower(double *L,int N,double fi,double fsval,double *dx) {
	int i,j,ct;
	double temp;
	
	if (fabs(fi) > fsval) {
		temp = fabs(fi);
	} else {
		temp = fsval;
	}
	
	temp = sqrt(temp);
	
	for(i = 0; i < N;++i) {
		ct = i *N;
		L[ct+i] = temp * dx[i];
		for(j = 0; j < i;++j) {
			L[ct+j] = 0;
		}
	}
	
}

void bfgs_factored(double *H,int N,double eps,double *xi,double *xf,double *jac,double *jacf) {
	int i,j,supd,ct;
	double sn,yn,fd,yt,jacm,alpha,temp3;
	double *sk,*yk,*temp,*temp2,*t,*u;
	
	sk = (double*) malloc(sizeof(double) *N);
	yk = (double*) malloc(sizeof(double) *N);
	temp = (double*) malloc(sizeof(double) *1);
	temp2 = (double*) malloc(sizeof(double) *1);
	t = (double*) malloc(sizeof(double) *N);
	u = (double*) malloc(sizeof(double) *N);
	
	msub(xf,xi,sk,1,N);
	msub(jacf,jac,yk,1,N);
	
	sn = l2norm(sk,N);
	yn = l2norm(yk,N);
	fd = sqrt(eps);
	
	mmult(yk,sk,temp,1,N,1);
	
	if (temp[0] >= fd*sn*yn) {
		for(i = 0; i < N;++i) {
			t[i] = 0.0;
			for(j = i; j < N;++j) {
				ct = j * N;
				t[i] += H[ct + i] * sk[j];
			}
		}
		mmult(t,t,temp2,1,N,1);
		alpha = sqrt(temp[0]/temp2[0]);
		supd = 1;
		for(i = 0; i < N;++i) {
			temp3 = 0.0;
			ct = i * N;
			for(j = 0; j < i+1;++j) {
				temp3 += H[ct + j] * t[j];
			}
			yt = fabs(yk[i] - temp3);
			if (jac[i] > jacf[i]) {
				jacm = jac[i];
			} else {
				jacm = jacf[i];
			}
			if (yt >= fd * jacm) {
				supd = 0;
			}
			u[i] = yk[i] - alpha*temp3;
		}
		if (supd == 0) {
			temp3 = 1.0 / sqrt(temp[0] * temp2[0]);
			for(i = 0; i < N; ++i) {
				t[i] *= temp3;
			}
			for(i = 1; i < N; ++i) {
				ct = i *N;
				for(j = 0; j < i;++j) {
					H[j*N+i] = H[ct+j];
				}
			}
			qrupdate(H,N,t,u);
			for(i = 1; i < N; ++i) {
				ct = i *N;
				for(j = 0; j < i;++j) {
					H[ct+j] = H[j*N+i];
				}
			}
		} 
	} 
	
	free(sk);
	free(yk);
	free(temp);
	free(temp2);
	free(t);
	free(u);
}

int bfgs_min(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double fsval,double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double stol,double *xf)  {
	int rcode,gfdcode;
	int i,siter,retval;
	double dt1,dt2;
	double fx,num,den,stop0,fxf,eps2;
	double *jac,*hess,*scheck,*xc,*L,*step,*jacf;
	
	jac = (double*) malloc(sizeof(double) *N);
	scheck = (double*) malloc(sizeof(double) *N);
	xc = (double*) malloc(sizeof(double) *N);
	step = (double*) malloc(sizeof(double) *N);
	hess = (double*) malloc(sizeof(double) *N * N);
	L = (double*) malloc(sizeof(double) *N * N);
	jacf = (double*) malloc(sizeof(double) *N);
	
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
	 15 -  Failure as Inf/Nan Values encountered 
	 * 
	 */ 
	
	rcode = 0;
	*niter = 0;
	siter = MAXITER;
	eps2 = sqrt(eps);
	gfdcode = 0;

	//set values
	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	fx = FUNCPT_EVAL(funcpt, xi, N);
	if (fx >= DBL_MAX || fx <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		rcode = 15;
	}
	if (fx != fx) {
		printf("Program Exiting as the function returns NaN");
		rcode = 15;
	}

	gfdcode = grad_fd(funcpt,funcgrad,xi,N,dx,eps2,jac);
	if (gfdcode == 15) {
		rcode = 15;
	}


	//maxstep = 1000.0; // Needs to be set at a much higher value proportional to l2 norm of dx

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
	
	//Check Stop0
	if (fabs(fx) > fabs(fsval)) {
			den = fabs(fx);
	} else {
			den = fabs(fsval);
	}
	for(i = 0; i < N;++i) {
		if (fabs(xi[i]) > 1.0 / fabs(dx[i])) {
			num = fabs(xi[i]);
		} else {
			num = 1.0 / fabs(dx[i]);
		}
		scheck[i] = fabs(jac[i]) * num / den;
	}
	
	stop0 = array_max_abs(scheck,N);
	
	if (stop0 <= gtol * 1e-03) {
		rcode = 1;
		for(i = 0; i < N;++i) {
			xf[i] = xi[i];
		}
	}
	
	//hessian_fd(funcpt,xi,N,dx,hess);
	inithess_lower(L,N,fx,fsval,dx);
	
	for(i = 0; i < N;++i) {
		xc[i] = xi[i];
	}
	
	while (rcode == 0 && *niter < siter) {
		*niter = *niter + 1;
		scale(jac,1,N,-1.0);
		
		linsolve_lower(L,N,jac,step);
		
		scale(jac,1,N,-1.0);

		retval = lnsrch(funcpt,xc,jac,step,N,dx,maxstep,stol,xf);

		fxf = FUNCPT_EVAL(funcpt, xf, N);
		if (fxf >= DBL_MAX || fxf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			rcode = 15;
			break;
		}
		if (fxf != fxf) {
			printf("Program Exiting as the function returns NaN");
			rcode = 15;
			break;
		}

		gfdcode = grad_fd(funcpt,funcgrad,xf,N,dx,eps2,jacf);
		if (gfdcode == 15) {
			rcode = 15;
			break;
		}
		rcode = stopcheck(fxf,N,xc,xf,jacf,dx,fsval,gtol,stol,retval);
		//hessian_fd(funcpt,xf,N,dx,hess);
		//bfgs_naive(hess,N,xc,xf,jac,jacf);
		bfgs_factored(L,N,eps,xc,xf,jac,jacf);
		for(i = 0; i < N;++i) {
			xc[i] = xf[i];
			jac[i] = jacf[i];
		}
	}
	
	if (rcode == 0 && *niter >= siter) {
		rcode = 4;
	}
	
	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	
	free(jac);
	free(hess);
	free(scheck);
	free(xc);
	free(L);
	free(step);
	free(jacf);
	return rcode;
}

int bfgs_min2(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, int m, double *dx, double fsval, double maxstep, int MAXITER, int *niter,
	double eps, double gtol, double ftol, double xtol, double *xf)  {
	int rcode, gfdcode;
	int i, siter, retval;
	double dt1, dt2;
	double fx, num, den, stop0, fxf, eps2,fo,alpha;
	double *jac, *hess, *scheck, *xc, *L, *step, *jacf;

	jac = (double*)malloc(sizeof(double)*N);
	scheck = (double*)malloc(sizeof(double)*N);
	xc = (double*)malloc(sizeof(double)*N);
	step = (double*)malloc(sizeof(double)*N);
	hess = (double*)malloc(sizeof(double)*N * N);
	L = (double*)malloc(sizeof(double)*N * N);
	jacf = (double*)malloc(sizeof(double)*N);

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
	* 15 -Failure as Inf/Nan Values encountered
	*
	*/

	rcode = 0;
	*niter = 0;
	siter = MAXITER;
	eps2 = sqrt(eps);

	alpha = 1.0;
	gfdcode = 0;

	//set values
	for (i = 0; i < N; ++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	fx = FUNCPT_EVAL(funcpt, xi, N);
	if (fx >= DBL_MAX || fx <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		rcode = 15;
	}
	if (fx != fx) {
		printf("Program Exiting as the function returns NaN");
		rcode = 15;
	}

	fo = fx;

	gfdcode = grad_fd(funcpt, funcgrad, xi, N, dx, eps2, jac);
	if (gfdcode == 15) {
		rcode = 15;
	}


	//maxstep = 1000.0; // Needs to be set at a much higher value proportional to l2 norm of dx

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

	//Check Stop0
	if (fabs(fx) > fabs(fsval)) {
		den = fabs(fx);
	}
	else {
		den = fabs(fsval);
	}
	for (i = 0; i < N; ++i) {
		if (fabs(xi[i]) > 1.0 / fabs(dx[i])) {
			num = fabs(xi[i]);
		}
		else {
			num = 1.0 / fabs(dx[i]);
		}
		scheck[i] = fabs(jac[i]) * num / den;
	}

	stop0 = array_max_abs(scheck, N);

	if (stop0 <= gtol * 1e-03) {
		rcode = 1;
		for (i = 0; i < N; ++i) {
			xf[i] = xi[i];
		}
	}

	//hessian_fd(funcpt,xi,N,dx,hess);
	inithess_lower(L, N, fx, fsval, dx);

	for (i = 0; i < N; ++i) {
		xc[i] = xi[i];
	}
	fxf = fx;

	while (rcode == 0 && *niter < siter) {
		*niter = *niter + 1;
		scale(jac, 1, N, -1.0);

		linsolve_lower(L, N, jac, step);

		scale(jac, 1, N, -1.0);

		for (i = 0; i < N; ++i) {
			jacf[i] = jac[i];
		}

		//retval = lnsrch(funcpt, xc, jac, step, N, dx, maxstep, stol, xf);
		retval = lnsrchmt(funcpt,funcgrad, xc, &fxf, jac, &alpha, step, N, dx, maxstep,MAXITER,eps2,ftol, gtol, xtol, xf);

		//rcode = stopcheck(fxf, N, xc, xf, jacf, dx, fsval, gtol, stol, retval);
		rcode = stopcheck2_mt(fxf, N, fo, jac, dx, eps, gtol, ftol, retval);
		fo = fxf;
		//hessian_fd(funcpt,xf,N,dx,hess);
		//bfgs_naive(hess,N,xc,xf,jac,jacf);
		bfgs_factored(L, N, eps, xc, xf, jacf, jac);
		for (i = 0; i < N; ++i) {
			xc[i] = xf[i];
		}
	}

	if (rcode == 0 && *niter >= siter) {
		rcode = 4;
	}

	for (i = 0; i < N; ++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}

	free(jac);
	free(hess);
	free(scheck);
	free(xc);
	free(L);
	free(step);
	free(jacf);
	return rcode;
}

void bfgs_rec(double *H,int N,int iter,int m,double *jac,double *sk,double *yk,double *r) {
	int i,ptr,k,bound;
	double *q,*temp,*rho,*alpha,*temp1,*temp2,*qs,*beta;


	q = (double*) malloc(sizeof(double) *N);
	temp = (double*) malloc(sizeof(double) *1);
	temp1 = (double*) malloc(sizeof(double) *N);
	temp2 = (double*) malloc(sizeof(double) *N);
	alpha = (double*) malloc(sizeof(double) * m);
	beta = (double*) malloc(sizeof(double) * m);
	rho = (double*) malloc(sizeof(double) * m);
	qs = (double*) malloc(sizeof(double) * m);

	k = 0; ptr = 0; i = 0; bound = 0;
	
	for(i = 0;i < N;++i) {
		q[i] = jac[i];
		temp1[i] = 0.0;
		temp2[i] = 0.0;
	}

	if (iter <= m) {
		bound = iter;
	} else {
		bound = m;
	}

	for (i = 0; i < m;++i) {
		alpha[i] = beta[i] = 0.0;
	}

	for (i = 0; i < m;++i) {
		rho[i] = qs[i] = 0.0;
	}

	//mmult(yk,sk,temp,1,N,1);

	for (i = bound - 1; i >= 0; --i) {
		ptr = i * N;
		for(k = 0; k < N;++k) {
			temp1[k] = yk[ptr+k];
			temp2[k] = sk[ptr+k];
		}
		mmult(temp1,temp2,temp,1,N,1);
		rho[i] = 1.0/temp[0];
		mmult(temp2,q,alpha+i,1,N,1);
		alpha[i] *= rho[i];
		for(k = 0; k < N; ++k) {
			q[k] = q[k] - alpha[i] * yk[ptr + k];
		}
	}
	mmult(H,q,r,N,N,1);
	for (i = 0; i < bound;++i) {
		ptr = i * N;
		for(k = 0; k < N;++k) {
			temp1[k] = yk[ptr+k];
			//temp2[k] = r[ptr+k];
		}
		mmult(temp1,r,beta+i,1,N,1);
		beta[i] *= rho[i];
		for (k = 0; k < N;++k) {
			r[k] += sk[ptr+k] * (alpha[i] - beta[i]);
		}
	}

	

	free(q);
	free(temp);
	free(rho);
	free(alpha);
	free(temp1);
	free(temp2);
	free(qs);
	free(beta);
}

void inithess_l(double *H, int N, int k, double *tsk,double *tyk, double *dx) {
	int i, j, ct;
	double temp;
	double *num, *den;

	num = (double*)malloc(sizeof(double)*1);
	den = (double*)malloc(sizeof(double)*1);

	mmult(tsk, tyk, num, 1, N, 1);
	mmult(tyk, tyk, den, 1, N, 1);

	if (k > 0) {
		temp = num[0] / den[0];
	}
	else {
		temp = 1.0;
	}

	for (i = 0; i < N; ++i) {
		ct = i *N;
		for (j = 0; j < N; ++j) {
			if (i == j) {
				H[ct + j] = temp * dx[i] * dx[i];
			}
			else {
				H[ct + j] = 0.0;
			}

		}
	}

	free(num);
	free(den);

}

int bfgs_l_min(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, int m, double *dx, double fsval,double maxstep, int MAXITER, int *niter,
		double eps,double gtol,double ftol,double xtol,double *xf)  {
	int rcode,gfdcode;
	int i,j,siter,retval;
	int ptr,iter;
	double dt1,dt2,alpha,fo;
	double fx,num,den,stop0,fxf,eps2;
	double *jac,*scheck,*xc,*H,*step,*jacf,*sk,*yk,*tsk,*tyk;

	jac = (double*) malloc(sizeof(double) *N);
	scheck = (double*) malloc(sizeof(double) *N);
	xc = (double*) malloc(sizeof(double) *N);
	step = (double*) malloc(sizeof(double) *N);
	H = (double*) malloc(sizeof(double) *N * N);
	jacf = (double*) malloc(sizeof(double) *N);
	sk = (double*) malloc(sizeof(double) * m * N);
	yk = (double*) malloc(sizeof(double) * m * N);
	tsk = (double*) malloc(sizeof(double) *N);
	tyk = (double*) malloc(sizeof(double) *N);

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

	rcode = 0;
	*niter = 0;
	siter = MAXITER;
	eps2 = sqrt(eps);

	alpha = 1.0;
	gfdcode = 0;
	
	//set values
	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	fx = FUNCPT_EVAL(funcpt, xi, N);
	if (fx >= DBL_MAX || fx <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		rcode = 15;
	}
	if (fx != fx) {
		printf("Program Exiting as the function returns NaN");
		rcode = 15;
	}
	fo = fx;

	gfdcode = grad_fd(funcpt,funcgrad,xi,N,dx,eps2,jac);
	if (gfdcode == 15) {
		rcode = 15;
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


	//Check Stop0
	if (fabs(fx) > fabs(fsval)) {
			den = fabs(fx);
	} else {
			den = fabs(fsval);
	}
	for(i = 0; i < N;++i) {
		if (fabs(xi[i]) > 1.0 / fabs(dx[i])) {
			num = fabs(xi[i]);
		} else {
			num = 1.0 / fabs(dx[i]);
		}
		scheck[i] = fabs(jac[i]) * num / den;
	}

	stop0 = array_max_abs(scheck,N);

	if (stop0 <= gtol * 1e-03) {
		rcode = 1;
		for(i = 0; i < N;++i) {
			xf[i] = xi[i];
		}
	}

	//hessian_fd(funcpt,xi,N,dx,hess);
	//inithess_lower(L,N,fx,fsval,dx);

	for(i = 0; i < N;++i) {
		xc[i] = xi[i];
	}
	fxf = fx;

	for (i = 0; i < N; ++i) {
		tsk[i] = 0.0;
		tyk[i] = 0.0;
	}


	while (rcode == 0 && *niter < siter) {
		iter = *niter;
		inithess_l(H,N,iter,tsk,tyk,dx);

		bfgs_rec(H,N,iter,m,jac,sk,yk,step);
		//mdisplay(step,1,N);

		scale(step,1,N,-1.0);

		for (i = 0; i < N; ++i) {
			jacf[i] = jac[i];
		}

		//linsolve_lower(L,N,jac,step);

		//scale(jac,1,N,-1.0);
		//retval = lnsrchmod(funcpt,xc,jac,step,N,dx,maxstep,stol,xf,jacf);
		retval = lnsrchmt(funcpt,funcgrad, xc, &fxf, jac, &alpha, step, N, dx, maxstep,MAXITER,eps2,ftol, gtol, xtol, xf);
		//retval = lnsrch(funcpt,xc,jac,step,N,dx,maxstep,stol,xf);

		//retval = swolfe(funcpt,xc,jac,step,N,dx,maxstep,stol,xf);
		//printf("%g %g \n",fo,fxf);

		//grad_fd(funcpt,xf,N,dx,jacf);
		//rcode = stopcheck_mt(fxf, N, xc, xf, jac, dx, fsval, gtol, ftol, retval);
		rcode = stopcheck2_mt(fxf,N,fo,jac,dx,eps,gtol,ftol,retval);
		//rcode = stopcheck3_mt(xc,xf,fxf,N,fo,jac,dx,eps,gtol,ftol,retval);
		//printf("\n CODE %d", rcode);
		fo = fxf;
		for (i = 0; i < N;++i) {
			tsk[i] = xf[i] - xc[i];
			tyk[i] = jac[i] - jacf[i];
		}

		if (*niter >= m) {
			for (i = 0; i < (m-1) * N; ++i) {
				sk[i] = sk[N + i];
				yk[i] = yk[N + i];
			}
			j = 0;

			for (i = (m - 1) * N; i < m * N; ++i) {
				sk[i] = tsk[j];
				yk[i] = tyk[j];
				j++;
			}

		} else {
			ptr = *niter * N;

			for (i = 0; i < N;++i) {
				sk[ptr + i] = tsk[i];
				yk[ptr + i] = tyk[i];
			}

		}

		
		for(i = 0; i < N;++i) {
			xc[i] = xf[i];
			step[i] = 0.0;
		}
		*niter = *niter + 1;

	}

	//mdisplay(xc, 1, N);

	if (rcode == 0 && *niter >= siter) {
		rcode = 4;
	}
	
	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	


	free(jac);
	free(scheck);
	free(xc);
	free(H);
	free(step);
	free(jacf);
	free(sk);
	free(yk);
	free(tsk);
	free(tyk);
	return rcode;
}


