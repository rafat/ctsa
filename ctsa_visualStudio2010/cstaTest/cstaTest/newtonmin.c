#include "newtonmin.h"

/*
 * newtonmin.c
 *
 *  Copyright (c) 2014, Rafat Hussain
 *	License : BSD 3-Clause
 *	See COPYRIGHT for more details
 */
 
 // Reference - Dennis, Schnabel, Numerical Methods for Unconstrained Optimization and Non-Linear Equations

double cholmod(double *A, int N, double *L, double eps,double maxinp) {
	/*
	 * Algorithm 5.5.2 Dennis, Schnabel, Numerical Methods for Unconstrained Optimization and 
	 * Non-Linear Equations
	 */ 
	int i,j,k,step,step2;
	double ls;
	double *U22;
	double maxadd,maxdiag,minl,minl2,minljj;
	
	U22 = (double*) malloc(sizeof(double) * N * N);
	
	minl = sqrt(sqrt(eps)) * maxinp;
	maxadd = 0.0;
	
	for(i = 0; i < N;++i) {
		U22[i] = A[i+i*N];
	}
	
	for(i = 0; i < N * N;++i) {
		L[i] = 0.0;
	}
	
	maxdiag = array_max_abs(U22,N);
	
	if ( maxinp == 0.0) {
		maxinp = sqrt(maxdiag);
	}
	
	minl2 = sqrt(eps) * maxinp;
	
	for(j = 0; j < N;++j) {
		step = j * N;
		ls = 0.0;
		for(i = 0; i < j;++i) {
			ls += L[step+i] * L[step+i];
		}
		
		L[step+j] = A[step+j] - ls;
		minljj = 0.0;
		for(i = j+1; i < N;++i) {
			ls = 0.0;
			step2 = i * N;
			for(k = 0;k < j;++k) {
				ls += L[step+k] * L[step2+k];
			}
			L[step2+j] = A[step+i] - ls;
			
			if (fabs(L[step2+j]) > minljj) {
				minljj = fabs(L[step2+j]);
			}
		}
		
		if (minljj/maxinp > minl) {
			minljj = minljj/maxinp;
		} else {
			minljj = minl;
		}
		
		if (L[step+j] > minljj*minljj) {
				L[step+j] = sqrt(L[step+j]);
		} else {
			if (minljj < minl2) {
				minljj = minl2;
			}
			
			if (maxadd < (minljj*minljj - L[step+j])) {
				maxadd = minljj*minljj - L[step+j];
			}
			
			L[step+j] = minljj;
		}
		
		for(i = j+1; i < N;++i) {
			L[i*N+j] /= L[step+j];
		}
		
		
	}
	
	free(U22);
	
	return maxadd;
}

double modelhess(double *A,int N,double *dx,double eps,double *L) {
	/*
	 * Algorithm 5.5.1 Dennis, Schnabel, Numerical Methods for Unconstrained Optimization and 
	 * Non-Linear Equations
	 */ 
	 double *U22;
	 double sqrteps,maxdiag,mindiag,maxposdiag,u;
	 double maxoffdiag,maxinp,maxadd;
	 double maxev,minev,offrow,sdd;
	 int step,i,j,k;
	 
	 sqrteps = sqrt(eps);
	 
	 U22 = (double*) malloc(sizeof(double) * N * N);
	 
	//scale
	
	for(i = 0;i < N;++i) {
		step = i*N;
		for(j = 0;j < N;++j) {
			A[step+j] /= (dx[i] * dx[j]);
		}
	}
	
	for(i = 0; i < N;++i) {
		U22[i] = A[i+i*N];
	}
	
	maxdiag = array_max(U22,N);
	mindiag = array_min(U22,N);
	
	maxposdiag = 0.0;
	
	if (maxdiag > maxposdiag) {
		maxposdiag = maxdiag;
	}
	
	for(i = 0; i < N;++i) {
		U22[i] = 0.0;
	}
	u = 0.0;
	if (mindiag <= sqrteps*maxposdiag) {
		u = 2 * (maxposdiag - mindiag) * sqrteps - mindiag;
		maxdiag += u;
	}
	k = 0;
	for (i = 0; i < N;++i) {
		for(j = 0;j < N;++j) {
			if (j > i) {
				U22[k] = A[i*N+j];
				k++;
			}
		}
	}
	
	maxoffdiag = array_max_abs(U22,k);
	
	if ( maxoffdiag*(1+2*sqrteps) > maxdiag) {
		u += (maxoffdiag - maxdiag) + 2*sqrteps*maxoffdiag;
		maxdiag = maxoffdiag*(1+2*sqrteps);
	}
	
	if (maxdiag == 0) {
		u = 1;
		maxdiag = 1;
	}
	
	if (u > 0) {
		for(i=0;i < N;++i) {
			A[i*N+i] += u;
		}
	}
	
	if (maxdiag > maxoffdiag / N) {
		maxinp = sqrt(maxdiag);
	} else {
		maxinp = sqrt(maxoffdiag / N);
	}
	
	maxadd = cholmod(A,N,L,eps,maxinp);
	
	if (maxadd > 0) {
		maxev = minev = A[0];
		for(i = 0; i < N;++i) {
			offrow = 0.0;
			step = i*N;
			
			for(j = 0; j < i;++j) {
				offrow += fabs(A[step+j]);
			}
			
			for(j = i+1; j < N;++j) {
				offrow += fabs(A[step+j]);
			}
			
			if (maxev < A[step+i] + offrow) {
				maxev = A[step+i] + offrow;
			}
			if (minev > A[step+i] - offrow) {
				minev = A[step+i] - offrow;
			}
			
		}
		sdd = (maxev - minev) * sqrteps - minev;
		if (sdd < 0) {
			sdd = 0;
		}
		if (maxadd > sdd) {
			u = sdd;
		} else {
			u = maxadd;
		}
		
		for(i = 0; i < N;++i) {
			A[i*N+i] += u; 
		}
	}
	
	maxadd = cholmod(A,N,L,eps,0.0);
	
	//unscale
	
	for(i = 0;i < N;++i) {
		step = i*N;
		for(j = 0;j < N;++j) {
			A[step+j] *= (dx[i] * dx[j]);
		}
	}
	
	for(i = 0;i < N;++i) {
		step = i*N;
		for(j = 0;j < i+1;++j) {
			L[step+j] *= dx[i];
		}
	}
	
	free(U22);
	
	return maxadd;
	 
}

void linsolve_lower(double *L,int N,double *b,double *x) {
	int i,j,c1,l;
	double *y,*A;
	double sum;
	
	y = (double*) malloc(sizeof(double) *N);
	A = (double*) malloc(sizeof(double) *N *N);
	
	for(i = 0; i < N;++i) {
		y[i] = 0.;
		x[i] = 0.;
		if ( L[i*N + i] == 0.) {
			printf("The Matrix system does not have a unique solution");
			exit(1);
		}
		//printf("\n B %d",ipiv[i]);
	}
	
	// Forward Substitution
	
	y[0] = b[0]/L[0];
	for(i = 1; i < N; ++i) {
		sum = 0.;
		c1 = i*N;
		for(j = 0; j < i; ++j) {
			sum += y[j] * L[c1 + j];
		}
		y[i] = (b[i] - sum)/L[c1+i];
	}
	
	mtranspose(L,N,N,A);
	
	//Back Substitution
	
	x[N - 1] = y[N - 1]/A[N * N - 1];
	
	for (i = N - 2; i >= 0; i--) {
		sum = 0.;
		c1 = i*(N+1);
		l=0;
		for(j = i+1; j < N;j++) {
			l++;
			sum += A[c1 + l] * x[j];
		}
		x[i] = (y[i] - sum) / A[c1];
	}
	
	free(y);
	free(A);
}



int hessian_fd(custom_function *funcpt,double *x,int N,double *dx,double eps,double *f) {
	int i,j,retval;
	double stepi,stepj,fd,stepmax;
	double *xi,*xj,*xij;
	
	retval = 0;

	fd = pow((double) eps,1.0/3.0);
	xi = (double*) malloc(sizeof(double) *N);
	xj = (double*) malloc(sizeof(double) *N);
	xij = (double*) malloc(sizeof(double) *N);
	
	for(i = 0; i < N;++i) {
		if (fabs(x[i]) >= 1.0 / fabs(dx[i])) {
			stepmax = x[i];
		} else {
			stepmax = signx(x[i]) * 1.0 / fabs(dx[i]);
		}
		stepi = fd * stepmax;
		for(j = 0; j < N;++j) {
			xi[j] = x[j];
			xj[j]= x[j];
			xij[j] = x[j];
		}
		xi[i] += stepi;
		xij[i] += stepi;
		for(j = 0; j < N;++j) {
			if (fabs(x[j]) >= 1.0 / fabs(dx[j])) {
				stepmax = x[j];
			} else {
				stepmax = signx(x[j]) * 1.0 / fabs(dx[j]);
			}
			stepj = fd * stepmax;
			xj[j] += stepj;
			xij[j] += stepj;
			//printf("stepmax %g %g \n",stepi,stepj);

			f[i*N+j] = ((FUNCPT_EVAL(funcpt,xij,N) - FUNCPT_EVAL(funcpt,xi,N)) - (FUNCPT_EVAL(funcpt,xj,N) - FUNCPT_EVAL(funcpt,x,N)))/(stepi * stepj);
			if (f[i*N+j] >= DBL_MAX || f[i*N+j] <= -DBL_MAX) {
				printf("Program Exiting as the function value exceeds the maximum double value");
				free(xi);
				free(xj);
				free(xij);
				return 15;
			}
			if (f[i*N+j] != f[i*N+j]) {
				printf("Program Exiting as the function returns NaN");
				free(xi);
				free(xj);
				free(xij);
				return 15;
			}
			xj[j] -= stepj;
			xij[j] -= stepj;
		}

	}
	
	free(xi);
	free(xj);
	free(xij);
	return retval;
	
}

int hessian_fd2(custom_function *funcpt,double *x,int N,double *dx,double eps,double *f) {
	int i,j,step,retval;
	double fd,temp,ft,fc,temp2,stepmax;
	double *stepsize,*f2;
	/*
	Returns Upper Triangular matrix only. If you want the full matrix then either use hessian_fd or copy the upper triangular value to lower triangle.
	*/
	fd = pow((double) eps,1.0/3.0);
	
	stepsize = (double*) malloc(sizeof(double) *N);
	f2 = (double*) malloc(sizeof(double) *N);
	retval = 0;
	fc = FUNCPT_EVAL(funcpt,x,N);
	if (fc >= DBL_MAX || fc <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		free(stepsize);
		free(f2);
		return 15;
	}
	if (fc != fc) {
		printf("Program Exiting as the function returns NaN");
		free(stepsize);
		free(f2);
		return 15;
	}
	
	for(i = 0; i < N;++i) {
		if (fabs(x[i]) >= (1.0 / fabs(dx[i]))) {
			stepmax = x[i];
		} else {
			stepmax = signx(x[i]) * 1.0 / fabs(dx[i]);
		}
		stepsize[i] = fd * stepmax;
		temp = x[i];
		x[i] += stepsize[i];
		stepsize[i] = x[i] - temp;
		
		f2[i] = FUNCPT_EVAL(funcpt,x,N);
		if (f2[i] >= DBL_MAX || f2[i] <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			free(stepsize);
			free(f2);
			return 15;
		}
		if (f2[i] != f2[i]) {
			printf("Program Exiting as the function returns NaN");
			free(stepsize);
			free(f2);
			return 15;
		}
		x[i] = temp;
	}
	
	for(i = 0; i < N;++i) { 
		step = i *N;
		temp = x[i];
		x[i] += 2*stepsize[i];
		ft = FUNCPT_EVAL(funcpt, x, N);
		if (ft >= DBL_MAX || ft <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			free(stepsize);
			free(f2);
			return 15;
		}
		if (ft != ft) {
			printf("Program Exiting as the function returns NaN");
			free(stepsize);
			free(f2);
			return 15;
		}
		f[step+i] = (fc - f2[i] + ft - f2[i])/(stepsize[i] * stepsize[i]);
		x[i] = temp + stepsize[i];
		for(j = i+1; j < N;++j) {
			temp2 = x[j];
			x[j] += stepsize[j];
			ft = FUNCPT_EVAL(funcpt, x, N);
			if (ft >= DBL_MAX || ft <= -DBL_MAX) {
				printf("Program Exiting as the function value exceeds the maximum double value");
				free(stepsize);
				free(f2);
				return 15;
			}
			if (ft != ft) {
				printf("Program Exiting as the function returns NaN");
				free(stepsize);
				free(f2);
				return 15;
			}
			f[step+j] = (fc - f2[i] + ft - f2[j])/(stepsize[i] * stepsize[j]);
			x[j] = temp2;
		}
		x[i] = temp;
			
	}
	
	free(stepsize);
	free(f2);
	return retval;
}

void fdjac(custom_gradient *funcgrad, double *x, int N, double *jac, double *dx, double eps2, double *J) {
	int i,j;
	double stepsize,temp,stepmax;
	double *fj;

	fj = (double*) malloc(sizeof(double) *N);

	for(j = 0; j < N;++j) {

		if (fabs(x[j]) >= (1.0 / fabs(dx[j]))) {
			stepmax = x[j];
		} else {
			stepmax = signx(x[j]) * 1.0 / fabs(dx[j]);
		}

		stepsize = eps2 * stepmax;
		temp = x[j];
		x[j] += stepsize;
		stepsize = x[j] - temp;
		FUNCGRAD_EVAL(funcgrad, x, N, fj);

		for(i = 0;i < N;++i) {
			J[i*N+j] = (fj[i] - jac[i]) / stepsize;
		}
		x[j] = temp;

	}

	free(fj);
}


void hessian_fdg(custom_gradient *funcgrad, double *x, int N, double *jac, double *dx, double eps2, double *H) {
	int i,j,k;

	fdjac(funcgrad,x,N,jac,dx,eps2,H);

	for(i = 0; i < N;++i) {
		k = i * N;
		for(j = 0; j < N;++j) {
			H[k + j] = (H[k + j] + H[j*N + i]) / 2.0;
		}
	}
}

int hessian_opt(custom_function *funcpt, custom_gradient *funcgrad, double *x, int N, double *jac,
		double *dx,double eps,double eps2,double *H) {
	int retval;
	retval = 0;
	if (funcgrad == NULL) {
		//printf("HESSF \n");
		retval = hessian_fd(funcpt,x,N,dx,eps,H);
	} else {
		//printf("HESSG \n");
		hessian_fdg(funcgrad,x,N,jac,dx,eps2,H);
	}

	return retval;
}

int lnsrch(custom_function *funcpt,double *xi,double *jac,double *p,int N,double * dx,double maxstep,double stol,double *x) {
	int retval,i,iter,MAXITER;
	double alpha,lambda,lambdamin,funcf,funci,lambdaprev,lambdatemp,funcprev;
	double lambda2,lambdaprev2,ll,den,rell,nlen;
	double *slopei,*temp1,*temp2,*ab,*rcheck,*pl;
	
	slopei = (double*) malloc(sizeof(double) *1);
	temp1 = (double*) malloc(sizeof(double) *4);
	temp2 = (double*) malloc(sizeof(double) *2);
	ab = (double*) malloc(sizeof(double) *2);
	rcheck = (double*) malloc(sizeof(double) *N);
	pl = (double*) malloc(sizeof(double) *N);
	retval = 100;
	alpha = 1e-04;
	lambda = 1.0;
	nlen = 0.0;
	MAXITER = 1000;
	funcprev = 1.0; // funcprev and lambdaprev are initialized to suppress warnings
	lambdaprev = 1.0; // These values are not used as the program sets the value later on.
	for(i = 0; i < N;++i) {
		nlen += dx[i] * p[i] * dx[i] * p[i];
	}
	nlen = sqrt(nlen);
	iter = 0;
	if (nlen > maxstep) {
		scale(p,1,N,maxstep/nlen);
		nlen = maxstep;
	}
	
	mmult(jac,p,slopei,1,N,1);
	for(i = 0; i < N;++i) {
		if (fabs(xi[i]) > 1.0 /fabs(dx[i])) {
			den = fabs(xi[i]);
		} else {
			den = 1.0 /fabs(dx[i]);
		}
		rcheck[i] = p[i]/den;
	}
	
	rell = array_max_abs(rcheck,N);
	
	lambdamin = stol/rell;
	
	//mdisplay(p,1,N);
	funci = FUNCPT_EVAL(funcpt, xi, N);
	if (funci >= DBL_MAX || funci <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		free(slopei);
		free(temp1);
		free(temp2);
		free(ab);
		free(rcheck);
		free(pl);
		return 15;
	}
	if (funci != funci) {
		printf("Program Exiting as the function returns NaN");
		free(slopei);
		free(temp1);
		free(temp2);
		free(ab);
		free(rcheck);
		free(pl);
		return 15;
	}
	while (retval > 1 && iter < MAXITER) {
		iter++;
		for(i = 0; i < N;++i) {
			pl[i] = p[i] * lambda;
		}
		madd(xi,pl,x,1,N);
		funcf = FUNCPT_EVAL(funcpt, x, N);
		if (funcf >= DBL_MAX || funcf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			free(slopei);
			free(temp1);
			free(temp2);
			free(ab);
			free(rcheck);
			free(pl);
			return 15;
		}
		if (funcf != funcf) {
			printf("Program Exiting as the function returns NaN");
			free(slopei);
			free(temp1);
			free(temp2);
			free(ab);
			free(rcheck);
			free(pl);
			return 15;
		}

		//printf("lambda %g %g %g \n",funcf,funci,lambda);
		
		if (funcf <= funci + alpha *lambda *slopei[0]) {
			retval = 0;
		} else if (lambda < lambdamin) {
			retval = 1;
			for (i = 0; i < N;++i) {
				x[i] = xi[i]; // Check
			}
		} else {
			if (lambda == 1.0) {
				lambdatemp = - slopei[0] / (2.0 * (funcf - funci - slopei[0])); 
			} else {
				lambda2 = lambda * lambda;
				lambdaprev2 = lambdaprev * lambdaprev;
				ll = lambda - lambdaprev;
				temp1[0] = 1.0 / lambda2; temp1[1] = -1.0 /lambdaprev2;
				temp1[2] = - lambdaprev / lambda2; temp1[3] = lambda /lambdaprev2;
				temp2[0] = funcf - funci - lambda * slopei[0];
				temp2[1] = funcprev - funci - lambdaprev * slopei[0];
				mmult(temp1,temp2,ab,2,2,1);
				scale(ab,1,2,1.0/ll);
				if (ab[0] == 0.0) {
					lambdatemp = - slopei[0] / (2.0 * ab[1]);
				} else {
					lambdatemp = (-ab[1] + sqrt( ab[1] * ab[1] - 3.0 * ab[0] *slopei[0]))/ (3.0 * ab[0]);
				}
				
				if (lambdatemp > 0.5 * lambda) {
					lambdatemp = 0.5 * lambda;
				}
			}
			lambdaprev = lambda;
			funcprev = funcf;
			if (lambdatemp <= 0.1 * lambda) {
				lambda = 0.1 * lambda;
			} else {
				lambda = lambdatemp;
			}
		}
	
	}
	
	free(slopei);
	free(temp1);
	free(temp2);
	free(ab);
	free(rcheck);
	free(pl);
	return retval;
}

int lnsrchmod(custom_function *funcpt, custom_gradient *funcgrad, double *xi, double *jac, double *p, int N, double * dx, double maxstep,
		double eps2,double stol,double *x,double *jacf) {
	int retval,i,gfdcode;
	double alpha,lambda,lambdamin,funcf,funci,lambdaprev,lambdatemp,funcprev;
	double lambda2,lambdaprev2,ll,den,rell,nlen;
	double *slopei,*temp1,*temp2,*ab,*rcheck,*pl;
	double *slopen;
	double beta,lambdamax,lambdalo,lambdainc,lambdadiff,funclo,funchi;
	
	slopei = (double*) malloc(sizeof(double) *1);
	slopen = (double*) malloc(sizeof(double) *1);
	temp1 = (double*) malloc(sizeof(double) *4);
	temp2 = (double*) malloc(sizeof(double) *2);
	ab = (double*) malloc(sizeof(double) *2);
	rcheck = (double*) malloc(sizeof(double) *N);
	pl = (double*) malloc(sizeof(double) *N);
	retval = 100;
	alpha = 1e-04;
	lambda = 1.0;
	nlen = 0.0;
	beta = 0.4;
	gfdcode = 0;
	funcprev = 1.0; // funcprev and lambdaprev are initialized to suppress warnings
	lambdaprev = 1.0; // These values are not used as the program sets the value later on.
	
	for(i = 0; i < N;++i) {
		nlen += dx[i] * p[i] * dx[i] * p[i];
	}
	nlen = sqrt(nlen);
	
	if (nlen > maxstep) {
		scale(p,1,N,maxstep/nlen);
		nlen = maxstep;
	}
	//mdisplay(p,1,N);
	mmult(jac,p,slopei,1,N,1);
	for(i = 0; i < N;++i) {
		if (fabs(xi[i]) > 1.0 /fabs(dx[i])) {
			den = fabs(xi[i]);
		} else {
			den = 1.0 /fabs(dx[i]);
		}
		rcheck[i] = p[i]/den;
	}
	
	rell = array_max_abs(rcheck,N);
	
	lambdamin = stol/rell;
	
	funci = FUNCPT_EVAL(funcpt, xi, N);
	
	if (funci >= DBL_MAX || funci <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		free(slopei);
		free(temp1);
		free(temp2);
		free(ab);
		free(rcheck);
		free(pl);
		free(slopen);
		return 15;
	}
	if (funci != funci) {
		printf("Program Exiting as the function returns NaN");
		free(slopei);
		free(temp1);
		free(temp2);
		free(ab);
		free(rcheck);
		free(pl);
		free(slopen);
		return 15;
	}
	while (retval > 1) {
		for(i = 0; i < N;++i) {
			pl[i] = p[i] * lambda;
		}
		madd(xi,pl,x,1,N);
		funcf = FUNCPT_EVAL(funcpt, x, N);
		//printf("%g lmax %g %g \n",lambda,funcf,funci + alpha *lambda *slopei[0]);
		if (funcf >= DBL_MAX || funcf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			free(slopei);
			free(temp1);
			free(temp2);
			free(ab);
			free(rcheck);
			free(pl);
			free(slopen);
			return 15;
		}
		if (funcf != funcf) {
			printf("Program Exiting as the function returns NaN");
			free(slopei);
			free(temp1);
			free(temp2);
			free(ab);
			free(rcheck);
			free(pl);
			free(slopen);
			return 15;
		}
		if (funcf <= funci + alpha *lambda *slopei[0]) {
			gfdcode = grad_fd(funcpt,funcgrad,x,N,dx,eps2,jacf);
			if (gfdcode == 15) {
				return 15;
			}
			mmult(jacf,p,slopen,1,N,1);
			
			if(slopen[0] < beta * slopei[0]) {
					if (lambda == 1.0 && nlen < maxstep) {
						lambdamax = maxstep / nlen;
						
						while (funcf <= funci + alpha *lambda *slopei[0] && slopen[0] < beta * slopei[0] && lambda < lambdamax) {
							lambdaprev = lambda;
							funcprev = funcf;
							if ( 2*lambda < lambdamax) {
								lambda = 2*lambda;
							} else {
								lambda = lambdamax;
							}
							for(i = 0; i < N;++i) {
								pl[i] = p[i] * lambda;
							}
							madd(xi,pl,x,1,N);
							funcf = FUNCPT_EVAL(funcpt, x, N);
							if (funcf >= DBL_MAX || funcf <= -DBL_MAX) {
								printf("Program Exiting as the function value exceeds the maximum double value");
								free(slopei);
								free(temp1);
								free(temp2);
								free(ab);
								free(rcheck);
								free(pl);
								free(slopen);
								return 15;
							}
							if (funcf != funcf) {
								printf("Program Exiting as the function returns NaN");
								free(slopei);
								free(temp1);
								free(temp2);
								free(ab);
								free(rcheck);
								free(pl);
								free(slopen);
								return 15;
							}
							if (funcf <= funci + alpha *lambda *slopei[0]) {
								gfdcode = grad_fd(funcpt,funcgrad,x,N,dx,eps2,jacf);
								if (gfdcode == 15) {
									free(slopei);
									free(temp1);
									free(temp2);
									free(ab);
									free(rcheck);
									free(pl);
									free(slopen);
									return 15;
								}
								mmult(jacf,p,slopen,1,N,1);
							}
						}
						
					}
					if (lambda < 1.0 || (lambda > 1.0 && funcf > funci + alpha *lambda *slopei[0])) {
						
						if (lambda < lambdaprev) {
							lambdalo = lambda; 
						} else {
							lambdalo = lambdaprev; 
						}
						lambdadiff = fabs(lambdaprev - lambda);
						
						if (lambda < lambdaprev) {
							funclo = funcf;
							funchi = funcprev;
						} else {
							funchi = funcf;
							funclo = funcprev;
						}

						while ((slopen[0] < beta * slopei[0]) && lambdadiff >= lambdamin) {
							lambdainc = -slopen[0] * lambdadiff * lambdadiff / (2.0 * (funchi - (funclo + slopen[0] *lambdadiff)));
							if ( lambdainc < 0.2 * lambdadiff) {
								lambdainc = 0.2 * lambdadiff;
							}
							lambda = lambdalo + lambdainc;
							for(i = 0; i < N;++i) {
								pl[i] = p[i] * lambda;
							}
							madd(xi,pl,x,1,N);
							funcf = FUNCPT_EVAL(funcpt, x, N);
							if (funcf >= DBL_MAX || funcf <= -DBL_MAX) {
								printf("Program Exiting as the function value exceeds the maximum double value");
								free(slopei);
								free(temp1);
								free(temp2);
								free(ab);
								free(rcheck);
								free(pl);
								free(slopen);
								return 15;
							}
							if (funcf != funcf) {
								printf("Program Exiting as the function returns NaN");
								free(slopei);
								free(temp1);
								free(temp2);
								free(ab);
								free(rcheck);
								free(pl);
								free(slopen);
								return 15;
							}
							if (funcf > funci + alpha *lambda *slopei[0]) {
								lambdadiff = lambdainc;
								funchi = funcf;
							} else {
								gfdcode = grad_fd(funcpt,funcgrad,x,N,dx,eps2,jacf);
								if (gfdcode == 15) {
									free(slopei);
									free(temp1);
									free(temp2);
									free(ab);
									free(rcheck);
									free(pl);
									free(slopen);
									return 15;
								}
								mmult(jacf,p,slopen,1,N,1);
								if (slopen[0] < beta * slopei[0]) {
									lambdalo = lambda;
									lambdadiff -= lambdainc;
									funclo = funcf;
								}
							}
						}
						if (slopen[0] < beta * slopei[0]) {
							retval = -5;
							funcf = funclo;
							for(i = 0; i < N;++i) {
								pl[i] = p[i] * lambdalo;
							}
							madd(xi,pl,x,1,N);
						}
					}
						
			}
			
			if (retval != -5) {
				retval = 0;
			}
		} else if (lambda < lambdamin) {
			retval = 1;
			for (i = 0; i < N;++i) {
				x[i] = xi[i]; // Check
			}
		} else {
			if (lambda == 1.0) {
				lambdatemp = - slopei[0] / (2.0 * (funcf - funci - slopei[0])); 
			} else {
				lambda2 = lambda * lambda;
				lambdaprev2 = lambdaprev * lambdaprev;
				ll = lambda - lambdaprev;
				temp1[0] = 1.0 / lambda2; temp1[1] = -1.0 /lambdaprev2;
				temp1[2] = - lambdaprev / lambda2; temp1[3] = lambda /lambdaprev2;
				temp2[0] = funcf - funci - lambda * slopei[0];
				temp2[1] = funcprev - funci - lambdaprev * slopei[0];
				mmult(temp1,temp2,ab,2,2,1);
				scale(ab,1,2,1.0/ll);
				if (ab[0] == 0.0) {
					lambdatemp = - slopei[0] / (2.0 * ab[1]);
				} else {
					lambdatemp = (-ab[1] + sqrt( ab[1] * ab[1] - 3.0 * ab[0] *slopei[0]))/ (3.0 * ab[0]);
				}
				
				if (lambdatemp > 0.5 * lambda) {
					lambdatemp = 0.5 * lambda;
				}
			}
			lambdaprev = lambda;
			funcprev = funcf;
			if (lambdatemp <= 0.1 * lambda) {
				lambda = 0.1 * lambda;
			} else {
				lambda = lambdatemp;
			}
		}
	
	}
	
	if (retval == -5) {
		retval = 5;
	}
	
	free(slopei);
	free(temp1);
	free(temp2);
	free(ab);
	free(rcheck);
	free(pl);
	free(slopen);
	return retval;
}

int lnsrchcg(custom_function *funcpt, custom_gradient *funcgrad, double *xi, double *jac, double *p, int N, double * dx, double maxstep,
	double eps2,double stol,double *x,double *jacf) {
	int retval,i,gfdcode;
	double alpha,lambda,lambdamin,funcf,funci,lambdaprev,lambdatemp,funcprev;
	double lambda2,lambdaprev2,ll,den,rell,nlen;
	double *slopei,*temp1,*temp2,*ab,*rcheck,*pl;
	double *slopen;
	double beta,lambdamax,lambdalo,lambdainc,lambdadiff,funclo,funchi;
	
	slopei = (double*) malloc(sizeof(double) *1);
	slopen = (double*) malloc(sizeof(double) *1);
	temp1 = (double*) malloc(sizeof(double) *4);
	temp2 = (double*) malloc(sizeof(double) *2);
	ab = (double*) malloc(sizeof(double) *2);
	rcheck = (double*) malloc(sizeof(double) *N);
	pl = (double*) malloc(sizeof(double) *N);
	retval = 100;
	alpha = 1e-04;
	lambda = 1.0;
	nlen = 0.0;
	beta = 0.9;
	gfdcode = 0;
	funcprev = 1.0; // funcprev and lambdaprev are initialized to suppress warnings
	lambdaprev = 1.0; // These values are not used as the program sets the value later on.
	
	for(i = 0; i < N;++i) {
		nlen += dx[i] * p[i] * dx[i] * p[i];
	}
	nlen = sqrt(nlen);
	
	if (nlen > maxstep) {
		scale(p,1,N,maxstep/nlen);
		nlen = maxstep;
	}
	
	mmult(jac,p,slopei,1,N,1);
	for(i = 0; i < N;++i) {
		if (fabs(xi[i]) > 1.0 /fabs(dx[i])) {
			den = fabs(xi[i]);
		} else {
			den = 1.0 /fabs(dx[i]);
		}
		rcheck[i] = p[i]/den;
	}
	
	rell = array_max_abs(rcheck,N);
	
	lambdamin = stol/rell;
	
	funci = FUNCPT_EVAL(funcpt, xi, N);
	
	if (funci >= DBL_MAX || funci <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		return 15;
	}
	if (funci != funci) {
		printf("Program Exiting as the function returns NaN");
		return 15;
	}
	while (retval > 1) {
		for(i = 0; i < N;++i) {
			pl[i] = p[i] * lambda;
		}
		madd(xi,pl,x,1,N);
		funcf = FUNCPT_EVAL(funcpt, x, N);
		if (funcf >= DBL_MAX || funcf <= -DBL_MAX) {
			printf("Program Exiting as the function value exceeds the maximum double value");
			return 15;
		}
		if (funcf != funcf) {
			printf("Program Exiting as the function returns NaN");
			return 15;
		}
		if (funcf <= funci + alpha *lambda *slopei[0]) {
			gfdcode = grad_fd(funcpt,funcgrad,x,N,dx,eps2,jacf);
			if (gfdcode == 15) {
				return 15;
			}
			mmult(jacf,p,slopen,1,N,1);
			if(fabs(slopen[0]) >= - beta * slopei[0]) {
					if (lambda == 1.0 && nlen < maxstep) {
						lambdamax = maxstep / nlen;
						while (funcf <= funci + alpha *lambda *slopei[0] || fabs(slopen[0]) >= -beta * slopei[0] || lambda < lambdamax) {
							lambdaprev = lambda;
							funcprev = funcf;
							if ( 2*lambda < lambdamax) {
								lambda = 2*lambda;
							} else {
								lambda = lambdamax;
							}
							for(i = 0; i < N;++i) {
								pl[i] = p[i] * lambda;
							}
							madd(xi,pl,x,1,N);
							funcf = FUNCPT_EVAL(funcpt, x, N);
							if (funcf >= DBL_MAX || funcf <= -DBL_MAX) {
								printf("Program Exiting as the function value exceeds the maximum double value");
								return 15;
							}
							if (funcf != funcf) {
								printf("Program Exiting as the function returns NaN");
								return 15;
							}
							if (funcf <= funci + alpha *lambda *slopei[0]) {
								gfdcode = grad_fd(funcpt,funcgrad,x,N,dx,eps2,jacf);
								if (gfdcode == 15) {
									return 15;
								}
								mmult(jacf,p,slopen,1,N,1);
							}
						}
					}
					if (lambda < 1.0 || (lambda > 1.0 && funcf > funci + alpha *lambda *slopei[0])) {
						if (lambda < lambdaprev) {
							lambdalo = lambda; 
						} else {
							lambdalo = lambdaprev; 
						}
						lambdadiff = fabs(lambdaprev - lambda);
						
						if (lambda < lambdaprev) {
							funclo = funcf;
							funchi = funcprev;
						} else {
							funchi = funcf;
							funclo = funcprev;
						}
						while (fabs(slopen[0]) >= -beta * slopei[0] || lambdadiff > lambdamin) {
							lambdainc = -slopen[0] * lambdadiff * lambdadiff / (2.0 * (funchi - (funclo + slopen[0] *lambdadiff)));
							if ( lambdainc < 0.2 * lambdadiff) {
								lambdainc = 0.2 * lambdadiff;
							}
							lambda = lambdalo + lambdainc;
							for(i = 0; i < N;++i) {
								pl[i] = p[i] * lambda;
							}
							madd(xi,pl,x,1,N);
							funcf = FUNCPT_EVAL(funcpt, x, N);
							if (funcf >= DBL_MAX || funcf <= -DBL_MAX) {
								printf("Program Exiting as the function value exceeds the maximum double value");
								return 15;
							}
							if (funcf != funcf) {
								printf("Program Exiting as the function returns NaN");
								return 15;
							}
							if (funcf > funci + alpha *lambda *slopei[0]) {
								lambdadiff = lambdainc;
								funchi = funcf;
							} else {
								gfdcode = grad_fd(funcpt,funcgrad,x,N,dx,eps2,jacf);
								if (gfdcode == 15) {
									return 15;
								}
								mmult(jacf,p,slopen,1,N,1);
								if (fabs(slopen[0]) >= -beta * slopei[0]) {
									lambdalo = lambda;
									lambdadiff -= lambdainc;
									funclo = funcf;
								}
							}
						}
						if (fabs(slopen[0]) >= -beta * slopei[0]) {
							funcf = funclo;
							for(i = 0; i < N;++i) {
								pl[i] = p[i] * lambdalo;
							}
							madd(xi,pl,x,1,N);
						}
					}
						
			}
			
			retval = 0;
		} else if (lambda < lambdamin) {
			retval = 1;
			for (i = 0; i < N;++i) {
				x[i] = xi[i]; // Check
			}
		} else {
			if (lambda == 1.0) {
				lambdatemp = - slopei[0] / (2.0 * (funcf - funci - slopei[0])); 
			} else {
				
				lambda2 = lambda * lambda;
				lambdaprev2 = lambdaprev * lambdaprev;
				ll = lambda - lambdaprev;
				temp1[0] = 1.0 / lambda2; temp1[1] = -1.0 /lambdaprev2;
				temp1[2] = - lambdaprev / lambda2; temp1[3] = lambda /lambdaprev2;
				temp2[0] = funcf - funci - lambda * slopei[0];
				temp2[1] = funcprev - funci - lambdaprev * slopei[0];
				mmult(temp1,temp2,ab,2,2,1);
				scale(ab,1,2,1.0/ll);
				if (ab[0] == 0.0) {
					lambdatemp = - slopei[0] / (2.0 * ab[1]);
				} else {
					lambdatemp = (-ab[1] + sqrt( ab[1] * ab[1] - 3.0 * ab[0] *slopei[0]))/ (3.0 * ab[0]);
				}
				if (lambdatemp > 0.5 * lambda) {
					lambdatemp = 0.5 * lambda;
				}
			}
			lambdaprev = lambda;
			funcprev = funcf;
			if (lambdatemp <= 0.1 * lambda) {
				lambda = 0.1 * lambda;
			} else {
				lambda = lambdatemp;
			}
		}
	
	}
	
	free(slopei);
	free(temp1);
	free(temp2);
	free(ab);
	free(rcheck);
	free(pl);
	free(slopen);
	return retval;
}

int stopcheck(double fx,int N,double *xc,double *xf,double *jac,double *dx,double fsval,double gtol,double stol,int retval) {
	int rcode,i;
	double num,den;
	double stop0;
	double *scheck;
	
	rcode = 0;	
	if (retval == 1) {
		rcode = 3;
		return rcode;
	}
	if (retval == 5) {
		rcode = 5;
		return rcode;
	}
	if (retval == 15) {
		rcode = 15;
		return rcode;
	}
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

int newton_min_func(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double fsval, double maxstep, int MAXITER,
		int *niter,double eps,double gtol,double stol,double *xf) {
	int rcode,gfdcode,hdcode;
	int i,siter,retval;
	double dt1,dt2,eps2;
	double fx,num,den,stop0,fxf;
	double *jac,*hess,*scheck,*xc,*L,*step;
	
	jac = (double*) malloc(sizeof(double) *N);
	scheck = (double*) malloc(sizeof(double) *N);
	xc = (double*) malloc(sizeof(double) *N);
	step = (double*) malloc(sizeof(double) *N);
	hess = (double*) malloc(sizeof(double) *N * N);
	L = (double*) malloc(sizeof(double) *N * N);
	
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
	gfdcode = 0;
	hdcode = 0;
	//gtol = pow((double) FDVAL,1.0/3.0);
	
	//stol = gtol * gtol;
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
	
	//printf("dt1 dt2 %g \n", maxstep);
	//maxstep = 1e15;
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
	
	//hessian_fd(funcpt,xi,N,dx,eps2,hess);
	hdcode = hessian_opt(funcpt,funcgrad,xi,N,jac,dx,eps,eps2,hess);
	if (hdcode == 15) {
		rcode = 15;
	}
	
	for(i = 0; i < N;++i) {
		xc[i] = xi[i];
	}
	
	while (rcode == 0 && *niter < siter) {
		*niter = *niter + 1;
		modelhess(hess,N,dx,eps,L);

		scale(jac,1,N,-1.0);
		//mdisplay(hess,N,N);
		
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
		//printf("%d \n",iter);
		gfdcode = grad_fd(funcpt,funcgrad,xf,N,dx,eps2,jac);
		if (gfdcode == 15) {
			rcode = 15;
			break;
		}
		rcode = stopcheck(fxf,N,xc,xf,jac,dx,fsval,gtol,stol,retval);
		//hessian_fd(funcpt,xf,N,dx,eps,hess);
		hdcode = hessian_opt(funcpt,funcgrad,xf,N,jac,dx,eps,eps2,hess);
		if (hdcode == 15) {
			rcode = 15;
		}
		for(i = 0; i < N;++i) {
			xc[i] = xf[i];
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
	return rcode;
}


int trsrch(custom_function *funcpt,double *xi,double *jac,double *sN,int N,double * dx,double maxstep,
		int iter,double *L,double *hess,double stol,double *ioval,double eps,double *x) {
	int retval,i,j,method;		
	double alpha,beta;
	double nlen,t1,funcprev;
	double *temp,*step,*xprev;
	/*
	 * ioval[0] = delta; [trust region current radius]
	 * ioval[1] = lambda;
	 * ioval[2] = deltaprev; [trust region previous radius]
	 * ioval[3] = phi;
	 * ioval[4] = phi';
	 * ioval[5] = phi_init';
	 * ioval[5] and ioval[6] are "binary" flags.
	 */ 
	retval = 4;
	nlen = 0.0;
	ioval[7] = 1.0;
	method = 0;
	temp = (double*) malloc(sizeof(double) *N);
	step = (double*) malloc(sizeof(double) *N);
	xprev = (double*) malloc(sizeof(double) *N);

	funcprev = 0.0;
	
	for(i = 0; i < N;++i) {
		nlen += dx[i] * sN[i] * dx[i] * sN[i];
	}
	nlen = sqrt(nlen);
	if (iter == 1 || ioval[0] == -1) {
		ioval[1] = 0.;
		for (i = 0; i < N;++i) {
			temp[i] = jac[i] / dx[i];
		}
		alpha = l2norm(temp,N);
		alpha = alpha*alpha;
		
		beta = 0;
		for(i = 0; i < N;++i) {
			t1 = 0.;
			for(j = i; j < N;++j) {
				t1+= L[j*N+i]*jac[j]/(dx[j] * dx[j]);
			}
			beta += t1*t1;
		}
		ioval[0] = alpha * sqrt(alpha) / beta;
		
		if (ioval[0] > maxstep) {
			ioval[0] = maxstep;
		}
	}
	
	while (retval >= 2 && retval != 15) {
		trstep(jac,sN,N,dx,L,hess,nlen,ioval,eps,step);
		ioval[2] = ioval[0];
		retval = trupdate(funcpt,xi,jac,step,N,dx,maxstep,retval,L,hess,stol,method,ioval,xprev,&funcprev,x);
	}
	 
	free(temp); 
	free(step);
	free(xprev);
	return retval;

}

void trstep(double *jac,double *sN,int N,double * dx,double *L,double *hess,double nlen,double *ioval,double eps,
		double *step) {
	int i,c1,iter,j;
	double sum,tvec,l1,l2,maxadd,steplen;
	double hi,lo,lambdalo,lambdahi;
	double delta,lambda,deltaprev,phi,phid,phidinit;
	double *temp1,*temp2;
	
	hi = 1.5;lo = 0.75;
	delta = ioval[0];lambda = ioval[1]; deltaprev = ioval[2];
	phi = ioval[3]; phid = ioval[4] ; phidinit = ioval[5];
	temp1 = (double*) malloc(sizeof(double) *N);
	temp2 = (double*) malloc(sizeof(double) *N);
	
	if (nlen <= hi*delta) {
		ioval[6] = 1.0;
		for(i=0;i <N;++i) {
			step[i] = sN[i];
		}
		lambda = 0.;
		if (delta > nlen) {
			delta = nlen;
		}
		
	} else {
			ioval[6] = 0.0;
			if (lambda > 0.) {
				lambda -= (phi+deltaprev) * ((deltaprev - delta) + phi) /(delta * phid);
			}
			phi = nlen - delta;
			if (ioval[7] == 1.0) {
				ioval[7] = 0.0;
				for(i=0;i < N;++i) {
					temp1[i] = dx[i] * dx[i] * sN[i];
				}
				temp2[0] = temp1[0]/L[0];
				for(i = 1; i < N; ++i) {
					sum = 0.;
					c1 = i*N;
					for(j = 0; j < i; ++j) {
						sum += temp2[j] * L[c1 + j];
					}
					temp2[i] = (temp1[i] - sum)/L[c1+i];
				}
				tvec = l2norm(temp2,N);
				phidinit = - tvec*tvec/nlen;
			}
			lambdalo = - phi/phidinit;
			for(i=0;i < N;++i) {
				temp1[i] = jac[i]/dx[i];
			}
			tvec = l2norm(temp1,N);
			lambdahi = tvec / delta;
			iter = 0;
			while (iter == 0) {
				if (lambda < lambdalo || lambda > lambdahi) {
					l1 = sqrt(lambdalo * lambdahi);
					l2 = 0.001 *lambdahi;
					if (l1 > l2) {
						lambda = l1;
					} else {
						lambda = l2;
					}
				}
				
				for(i = 0; i < N;++i) {
					c1 = i*N;
					hess[c1+i] += lambda * dx[i] * dx[i];
				}
				maxadd = cholmod(hess,N,L,eps,0.0);
				scale(jac,1,N,-1.0);
				linsolve_lower(L,N,jac,step);
				scale(jac,1,N,-1.0);
				for(i = 0; i < N;++i) {
					c1 = i*N;
					hess[c1+i] -= lambda * dx[i] * dx[i];
				}
				for(i=0;i < N;++i) {
					temp1[i] = step[i] * dx[i];
				}
				steplen = l2norm(temp1,N);
				phi = steplen - delta;
				for(i=0;i < N;++i) {
					temp1[i] *= dx[i];
				}
				temp2[0] = temp1[0]/L[0];
				for(i = 1; i < N; ++i) {
					sum = 0.;
					c1 = i*N;
					for(j = 0; j < i; ++j) {
						sum += temp2[j] * L[c1 + j];
					}
					temp2[i] = (temp1[i] - sum)/L[c1+i];
				}
				tvec = l2norm(temp2,N);
				phid = - tvec*tvec/steplen;
				
				if ( (steplen >= lo * delta && steplen <= hi * delta ) || (lambdahi - lambdalo >= 0))  {
					iter = 1;
				} else {
					if (lambdalo <= lambda - (phi/phid)) {
						lambdalo = lambda - (phi/phid);
					}
					if (phi < 0) {
						lambdalo = lambda;
					}
					
					lambda -= (steplen *phi)/(delta*phid);
				}
				
			}
			
	}
	
	ioval[0] = delta;ioval[1] =lambda; ioval[2] = deltaprev;
	ioval[3] = phi; ioval[4] = phid ;ioval[5] = phidinit;
	
	free(temp1);
	free(temp2);
}
	
int trupdate(custom_function *funcpt,double *xi,double *jac,double *step,int N,double * dx,double maxstep,
		int retcode,double *L,double *hess,double stol,int method,double *ioval,double *xprev,double *funcprev,double *x) {
			
	int retval;	
	int i,j,c1;
	double delta,alpha,nlen,rell,den,deltatemp;
	double funci,funcf,df,dfp,temp;
	double *slopei,*rcheck;
	/*Return Codes
	 * 
	 * 0 - x Accepted.
	 * 1 - x unsatisfactory but (x-xi) is smaller than stopping criteria stol
	 * 2 - f(x) large. Reduce delta
	 * 3 - f(x) sufficiently small. Double delta
	 * 
	 */ 
	 
	alpha = 1.0e-04;
	delta = ioval[0];
	nlen = 0.0;
	 
	slopei = (double*) malloc(sizeof(double) *1);
	rcheck = (double*) malloc(sizeof(double) *N);
 
	for(i = 0; i < N;++i) {
		nlen += dx[i] * step[i] * dx[i] * step[i];
	}
	nlen = sqrt(nlen);
	
	
	funci = FUNCPT_EVAL(funcpt, xi, N);
	if (funci >= DBL_MAX || funci <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		free(rcheck);
		free(slopei);
		return 15;
	}
	if (funci != funci) {
		printf("Program Exiting as the function returns NaN");
		free(rcheck);
		free(slopei);
		return 15;
	}
	
	madd(xi,step,x,1,N);
	
	funcf = FUNCPT_EVAL(funcpt, x, N);
	if (funcf >= DBL_MAX || funcf <= -DBL_MAX) {
		printf("Program Exiting as the function value exceeds the maximum double value");
		free(rcheck);
		free(slopei);
		return 15;
	}
	if (funcf != funcf) {
		printf("Program Exiting as the function returns NaN");
		free(rcheck);
		free(slopei);
		return 15;
	}
	
	df = funcf - funci;
	mmult(jac,step,slopei,1,N,1);
	 
	if (retcode != 3) {
		*funcprev = 0.;
	} 
	
	if ((df > alpha *slopei[0] || funcf >= *funcprev) && retcode == 3) {
		for(i = 0;i < N;++i) {
			x[i] = xprev[i];
		}
		retval = 0;
		funcf = *funcprev;
		delta /= 2.0;
	} else if (df >= alpha *slopei[0]) {
		for(i = 0; i < N;++i) {
			if (fabs(xi[i]) > 1.0 /fabs(dx[i])) {
				den = fabs(xi[i]);
			} else {
				den = 1.0 /fabs(dx[i]);
			}
			rcheck[i] = step[i]/den;
		}
	
		rell = array_max_abs(rcheck,N);
		
		if (rell < stol) {
			retval = 1;
			for(i = 0; i < N;++i) {
				x[i] = xi[i];
			}
		} else {
			retval = 2;
			deltatemp = - slopei[0] * nlen /(2*(df - nlen));
			if (deltatemp < 0.1 * delta) {
				delta *= 0.1;
			} else if (deltatemp > 0.5 * delta) {
				delta *= 0.5;
			} else {
				delta = deltatemp;
			}
		}
		
	} else {
		dfp = slopei[0];
		if (method == 0) {
			for (i = 0;i < N;++i) {
				//Hook
				c1 = i * N;
				temp = 0.5 * hess[c1+i] * step[i] * step[i];
				for (j = i+1; j < N;++j) {
					temp += hess[c1+j] * step[i] * step[j];
				}
				dfp += temp;
			}
		} else if (method == 1) {
			//Double Dogleg
			for (i = 0; i < N;++i) {
				temp = 0.0;
				for(j = i; j < N;++j) {
					temp += L[j*N+i] * step[j];
				}
				dfp += temp*temp/2.0;
			}
		}
		
		if ((retcode != 2) && (fabs(dfp - df) <= fabs(df) || df <= slopei[0]) && (ioval[6] == 0.0) && (delta <= 0.99*maxstep)) {
			retval = 3;
			*funcprev = funcf;
			for(i = 0; i < N;++i) {
				xprev[i] = x[i];
			}
			if (maxstep > 2*delta) {
				delta = 2 * delta;
			} else {
				delta = maxstep;
			}
			
		} else {
			retval = 0;
			if (df >= 0.1 * dfp) {
				delta = delta/2;
			} else if (df <= 0.75 *dfp) {
				if (maxstep > 2*delta) {
				delta = 2 * delta;
				} else {
					delta = maxstep;
				}
				
			}
		}
	}
	
	
	ioval[0] = delta;
	free(rcheck);
	free(slopei); 
	return retval;
		
}

int newton_min_trust(custom_function *funcpt, custom_gradient *funcgrad, double *xi, int N, double *dx, double fsval, double delta,
		int method,int MAXITER,int *niter,double eps,double gtol,double stol,double *xf) {
	int rcode,iter,gfdcode,hdcode;
	int i,siter,retval,fdiff;
	double dt1,dt2,eps2,eps3;
	double fx,num,den,stop0,maxstep,fxf;
	double *jac,*hess,*scheck,*xc,*L,*step;
	double *ioval;
	
	jac = (double*) malloc(sizeof(double) *N);
	scheck = (double*) malloc(sizeof(double) *N);
	xc = (double*) malloc(sizeof(double) *N);
	step = (double*) malloc(sizeof(double) *N);
	hess = (double*) malloc(sizeof(double) *N * N);
	L = (double*) malloc(sizeof(double) *N * N);
	ioval = (double*) malloc(sizeof(double) *8);
	
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
	eps3 = pow(eps,(double) 1.0/3.0);
	eps2 = sqrt(eps);
	gfdcode = 0;
	hdcode = 0;
	//set values
	
	if (delta <= 0.0) {
		delta = -1.0;
	}
	ioval[0] = delta;
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
	fdiff = 1;
	//grad_cd(funcpt,funcgrad,xi,N,dx,eps3,jac);
	
	maxstep = 1000.0; // Needs to be set at a much higher value proportional to l2 norm of dx
	dt1 = dt2 = 0.0;
	for(i = 0; i < N;++i) {
		dt1 += dx[i] * dx[i];
		dt2 += dx[i] * xi[i] * dx[i] * xi[i];
	}

	dt1 = sqrt(dt1);
	dt2 = sqrt(dt2);
	
	if (dt1 > dt2) {
		maxstep *= dt1;
	} else {
		maxstep *= dt2;
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
	
	//hessian_fd(funcpt,xi,N,dx,eps,hess);
	hdcode = hessian_opt(funcpt,funcgrad,xi,N,jac,dx,eps,eps2,hess);
	if (hdcode == 15) {
		rcode = 15;
	}
	
	for(i = 0; i < N;++i) {
		xc[i] = xi[i];
	}
	
	while (rcode == 0 && *niter < siter) {
		*niter = *niter + 1;
		
		modelhess(hess,N,dx,eps,L);
		scale(jac,1,N,-1.0);
		
		linsolve_lower(L,N,jac,step);
		
		scale(jac,1,N,-1.0);
		iter = *niter;
		if (method == 0) {
			retval = trsrch(funcpt,xc,jac,step,N,dx,maxstep,iter,L,hess,stol,ioval,eps,xf);
		} else if (method == 1) {
			retval = trsrch_ddl(funcpt,xc,jac,step,N,dx,maxstep,iter,L,hess,stol,ioval,xf);
		}
		else {
			printf("The program accepts only two method values \n");
			printf("Method 0 : Hook step \n");
			printf("Method 1 : Double Dog Leg Step \n");
			exit(1);
		}
		
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

		gfdcode = grad_fd(funcpt,funcgrad,xf,N,dx,eps2,jac);
		if (gfdcode == 15) {
			rcode = 15;
			break;
		}
		//grad_cd(funcpt,funcgrad,xf,N,dx,eps3,jac);
		rcode = stopcheck(fxf,N,xc,xf,jac,dx,fsval,gtol,stol,retval);
		//hessian_fd(funcpt,xf,N,dx,eps,hess);
		hdcode = hessian_opt(funcpt,funcgrad,xf,N,jac,dx,eps,eps2,hess);
		if (hdcode == 15) {
			rcode = 15;
		}
		for(i = 0; i < N;++i) {
			xc[i] = xf[i];
		}
	}
	
	if (rcode == 0 && *niter >= siter) {
		rcode = 4;
	}
	
	for(i = 0; i < N;++i) {
		xi[i] *= dx[i];
		dx[i] = 1.0 / dx[i];
	}
	
	free(ioval);
	free(jac);
	free(hess);
	free(scheck);
	free(xc);
	free(L);
	free(step);
	return rcode;
}

int trsrch_ddl(custom_function *funcpt,double *xi,double *jac,double *sN,int N,double * dx,double maxstep,
		int iter,double *L,double *hess,double stol,double *ioval,double *x) {
	int retval,i,method;		
	double nlen,funcprev;
	double *temp,*step,*xprev,*ssd,*v;
	/*
	 * ioval[0] = delta; [trust region current radius]
	 * ioval[1] = cauchy step length;
	 * ioval[2]-ioval[5] dummy constants;
	 * ioval[6] and ioval[7] are "binary" flags.
	 */ 
	retval = 4;
	nlen = 0.0;
	ioval[7] = 1.0;
	method = 1;
	temp = (double*) malloc(sizeof(double) *N);
	step = (double*) malloc(sizeof(double) *N);
	xprev = (double*) malloc(sizeof(double) *N);
	ssd = (double*) malloc(sizeof(double) *N);
	v = (double*) malloc(sizeof(double) *N);
	
	for(i = 0; i < N;++i) {
		nlen += dx[i] * sN[i] * dx[i] * sN[i];
	}
	nlen = sqrt(nlen);
	funcprev = 0.0;
	
	while (retval >= 2 && retval != 15) {
		trstep_ddl(jac,sN,N,dx,maxstep,L,hess,nlen,ioval,ssd,v,step);
		retval = trupdate(funcpt,xi,jac,step,N,dx,maxstep,retval,L,hess,stol,method,ioval,xprev,&funcprev,x);
	}
	 
	free(temp); 
	free(step);
	free(xprev);
	free(ssd);
	free(v);
	return retval;

}
		
void trstep_ddl(double *jac,double *sN,int N,double * dx,double maxstep,double *L,double *hess,double nlen,double *ioval,
		double *ssd,double *v,double *step) {
	int i,j;
	double tvec,temp;
	double alpha,beta;
	double delta,clen,neta,lambda;
	double *temp1,*slstep;
	double *t1,*t2;
	
	
	delta = ioval[0];clen = ioval[1]; neta = ioval[2];
	
	temp1 = (double*) malloc(sizeof(double) *N);
	slstep = (double*) malloc(sizeof(double) *1);
	t1 = (double*) malloc(sizeof(double) *1);
	t2 = (double*) malloc(sizeof(double) *1);
	
	if (nlen <= delta) {
		ioval[6] = 1.0;
		for(i=0;i <N;++i) {
			step[i] = sN[i];
		}
		delta = nlen;
		
	} else {
		ioval[6] = 0.0;
		if (ioval[7] == 1.0) {
			ioval[7] = 0.0;
			for(i=0;i < N;++i) {
				temp1[i] = jac[i] / dx[i];
			}
			
			tvec = l2norm(temp1,N);
			alpha = tvec*tvec;
			beta = 0.0; 
			
			for (i = 0; i < N;++i) {
				temp = 0.0;
				for(j = i; j < N;++j) {
					temp += L[j*N + i] * jac[j] /(dx[j] * dx[j]);
				}
				beta += temp*temp;
			}
			tvec = - alpha / beta;
			for (i = 0; i < N;++i) {
				ssd[i] = tvec * jac[i] / dx[i];
			}
			
			clen = alpha * sqrt(alpha) / beta;
			mmult(jac,sN,slstep,1,N,1);
			neta = 0.2 + ((0.8 * alpha * alpha) /(beta *fabs(slstep[0])));
			
			for (i = 0; i < N;++i) {
				v[i] = neta * dx[i] *sN[i] - ssd[i];
			}
			
			if (delta == -1.0) {
				if (clen < maxstep) {
					delta = clen;
				} else {
					delta = maxstep;
				}
			}
			
		}
		
		if (neta * nlen <= delta) {
			tvec = delta / nlen;
			for(i = 0; i < N;++i) {
				step[i] = tvec * sN[i];
			}
		} else if (clen >= delta) {
			tvec = delta /clen;
			for (i = 0; i < N;++i) {
				step[i] = tvec * ssd[i] / dx[i];
			}
		} else {
			mmult(v,ssd,t1,1,N,1);
			mmult(v,v,t2,1,N,1);
			lambda = (-t1[0] + sqrt(t1[0]*t1[0] - t2[0]*(clen*clen - delta*delta)))/t2[0];
			
			for (i = 0; i < N; ++i) {
				step[i] = (ssd[i] / dx[i]) + (lambda * v[i] / dx[i] );
			}
		}
		
			
	}
	
	ioval[0] = delta; ioval[1] = clen; ioval[2] = neta;
	
	free(temp1);
	free(t1);
	free(t2);
	free(slstep);
}		
