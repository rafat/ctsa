/*
 * initest.c
 *
 *  Created on: Jul 13, 2014
 *      Author: Rafat Hussain
 */

#include "initest.h"

void innalg(double *vec, int l, int N,double *v, double *theta) {
	int i,j,lag,ct,ct2,index1,index2;
	double temp,temp2;
	double *acov;
	// IMPORTANT - The program does not store the matrix as laid out
	// in Brockwell and Davis

	// theta is a (N+1)X(N+1) matrix so choose the value of N accordingly. The recommended
	// value is N = 17 which will make theta a 18*18 double vector so allocate memory
	// accordingly.

	// theta11 is accessed at theta[1*(N+1) + 1]

	// with thetaxy accessed at theta[x*(N+1)+y]

	// The coefficients of Moving Average MA(q) are stored at
	// theta[N*(N+1)+1+i] where i = 0,1,2,...,q-1. Multiply these entries by -1.0
	// to get the MA(q) coefficients.

	// Output v is a 1X(N+1) vector with v[0] = autocovariance at lag 0.

	// see ma_inn function below

	lag = N+1;

	acov = (double*) malloc(sizeof(double) * lag);
	autocovar(vec,l,acov,lag);
	//mdisplay(acov,1,lag);
	v[0] = acov[0];

	for(i = 0; i < lag;++i) {
		index1 = i * lag;
		for(j = 0; j < lag;++j) {
			theta[index1+j] = 0.0;
		}
	}

	for(i = 1; i < lag;++i) {
		index1 = i * lag;
		for(j = 0; j < i;++j) {
			index2 = j * lag;
			temp = 0.0;
			for(ct = 0; ct < j;++ct) {
				temp += (theta[index2 + j - ct] * theta[index1 + i - ct]* v[ct]);
			}
			theta[index1 + i - j] = (acov[i-j]-temp)*1.0/v[j];
		}
		temp2 = 0.0;

		for(ct2 = 0; ct2 < i;++ct2) {
			temp2 += (v[ct2]*theta[index1 + i - ct2]*theta[index1 + i - ct2]);
		}

		v[i]=v[0]-temp2;
	}


	free(acov);
}

void ma_inn(double *inp, int N,int q, double* theta, int lag) {

	// Moving Average Parameters for a MA(q) process using Innovations
	// Algorithm. lag is typically set at <= 17
	int i;
	double *v,*thetap;

	thetap = (double*) malloc(sizeof(double) * (lag+1) * (lag+1));
	v = (double*) malloc(sizeof(double) * (lag+1));
	innalg(inp,N,lag,v,thetap);
	for(i = 0; i < q; ++i) {
		theta[i] = -thetap[lag*(lag+1) + 1 + i];
	}

	free(v);
	free(thetap);
}

void burgalg(double *x, int N, int p, double *phi,double *var) {
	int i, j, N1, p1;
	double dk,u,temp1,temp2;
	double *fk, *bk, *an;

	/*
	This code is the C translation of Cedrick Collomb's C++ Implemnetation of Burg's Algorithm
	with some minor modifications
	A tutorial on linear prediction and the Burg method.
	http://www.emptyloop.com/technotes/
	*/

	N1 = N;
	N = N - 1;
	p1 = p + 1;
	fk = (double*)malloc(sizeof(double) * N1);
	bk = (double*)malloc(sizeof(double)* N1);
	an = (double*)malloc(sizeof(double)* p1);

	for (i = 0; i < N1; ++i) {
		fk[i] = bk[i] = x[i];
	}

	for (i = 0; i < p1; ++i) {
		an[i] = 0.0;
		//printf("%g ", an[i]);
	}

	an[0] = 1.0;

	dk = 0.0;

	for (i = 0; i < N1; ++i) {
		dk += ( 2.0 * fk[i] * fk[i]);
	}

	dk -= (fk[0] * fk[0] + bk[N] * bk[N]);

	for (i = 0; i < p; ++i) {
		u = 0.0;
		for (j = 0; j <= N - i - 1; ++j) {
			u += (fk[i + j + 1] * bk[j]);
		}

		u = -2.0 * u / dk;
		 
		for (j = 0; j <= (i + 1) / 2; ++j) {
			temp1 = an[j] + u * an[i + 1 - j];
			temp2 = an[i + 1 - j] + u * an[j];
			an[j] = temp1;
			an[i + 1 - j] = temp2;
		}

		for (j = 0; j <= N - i - 1; ++j) {
			temp1 = fk[i + j + 1] + u * bk[j];
			temp2 = bk[j] + u * fk[i + j + 1];
			fk[i + j + 1] = temp1;
			bk[j] = temp2;
		}
		*var = (1.0 - u*u) * dk;
		dk = *var - fk[i + 1] * fk[i + 1] - bk[N - i - 1] * bk[N - i - 1];
		*var = *var / (2 * (N - i));

	}

	for (i = 0; i < p; ++i) {
		phi[i] = -1.0 * an[i+1];
	}

	free(fk);
	free(bk);
	free(an);
}

void burgalg2(double *x, int N, int p, double *phi, double *var) {
	int i, j, N1, p1;
	double dk, u, temp1, temp2;
	double *fk, *bk, *an;

	N1 = N;
	N = N - 1;
	p1 = p + 1;

	fk = (double*)malloc(sizeof(double)* N1);
	bk = (double*)malloc(sizeof(double)* N1);
	an = (double*)malloc(sizeof(double)* p1);

	for (i = 0; i < N1; ++i) {
		fk[i] = x[i+1];
		bk[i] = x[i];
	}

	free(fk);
	free(bk);
	free(an);
}

void ywalg(double *x, int N, int p, double *phi) {
	int lag,k,j;
	double *acf,*z,*y,*temp1,*temp2;
	double a, b;

	lag = p + 1;
	acf = (double*)malloc(sizeof(double)* lag);
	y = (double*)malloc(sizeof(double)* lag);
	z = (double*)malloc(sizeof(double)* p);
	temp1 = (double*)malloc(sizeof(double)* p);
	temp2 = (double*)malloc(sizeof(double)* 1);

	autocorr(x, N, acf, lag);
	b = acf[0];
	a = y[0] = -acf[1];

	for (k = 0; k < p - 1; ++k) {
		b = (1 - a*a) * b;
		for (j = k; j >= 0; --j) {
			temp1[k-j] = acf[j+1];
		}
		mmult(temp1, y, temp2, 1, k + 1, 1);

		a = -(*temp2 + acf[k + 2]) / b;

		for (j = 0; j <= k; ++j) {
			z[j] = y[j] + a * y[k - j];
		}

		for (j = 0; j <= k; ++j) {
			y[j] = z[j];
		}
		y[k + 1] = a;

	}

	for (j = 0; j < p; ++j) {
		phi[j] = -1.0 * y[j];
	}


	free(acf); 
	free(z);
	free(y);
	free(temp1);
	free(temp2);
}

void ywalg2(double *x, int N, int p, double *phi,double *var) {
	int lag, k, j;
	double *acf, *z, *y, *temp1, *temp2;
	double a, b;
	// Same as ywalg except it also returns the variance value
	lag = p + 1;
	acf = (double*)malloc(sizeof(double)* lag);
	y = (double*)malloc(sizeof(double)* lag);
	z = (double*)malloc(sizeof(double)* p);
	temp1 = (double*)malloc(sizeof(double)* p);
	temp2 = (double*)malloc(sizeof(double)* 1);

	autocovar(x, N, acf, lag);
	*var = acf[0];
	acf[0] = 1.0;
	for (k = 1; k < lag; ++k) {
		acf[k] = acf[k] / *var;
	}
	b = acf[0];
	a = y[0] = -acf[1];

	for (k = 0; k < p - 1; ++k) {
		b = (1 - a*a) * b;
		for (j = k; j >= 0; --j) {
			temp1[k - j] = acf[j + 1];
		}
		mmult(temp1, y, temp2, 1, k + 1, 1);

		a = -(*temp2 + acf[k + 2]) / b;

		for (j = 0; j <= k; ++j) {
			z[j] = y[j] + a * y[k - j];
		}

		for (j = 0; j <= k; ++j) {
			y[j] = z[j];
		}
		y[k + 1] = a;

	}

	for (j = 0; j < p; ++j) {
		phi[j] = -1.0 * y[j];
	}
	mmult(acf + 1, phi, temp2, 1, p, 1);
	*var = *var * (1.0 - *temp2);

	free(acf);
	free(z);
	free(y);
	free(temp1);
	free(temp2);
}

void hralg(double *x, int N, int p,int q, double *phi,double *theta, double *var) {
	int m,i,t,pq,j,k;
	double wmean,sos;
	double *inp,*phim,*a,*z,*zt,*A,*b,*temp;
	int *ipiv;

	inp = (double*)malloc(sizeof(double)* N);

	if (p > q) {
		m = p;
	}
	else {
		m = q;
	}

	if (m < 10  && N > p+q+20) {
		m = 10;
	}
	else {
		m = m + 1;
	}

	pq = p + q;

	phim = (double*)malloc(sizeof(double)* m);
	a = (double*)malloc(sizeof(double)* N);
	z = (double*)malloc(sizeof(double)* (N - m - q) * pq);
	zt = (double*)malloc(sizeof(double)* (N - m - q) * pq);
	A = (double*)malloc(sizeof(double)* pq * pq);
	b = (double*)malloc(sizeof(double)* pq);
	temp = (double*)malloc(sizeof(double)* pq);
	ipiv = (int*)malloc(sizeof(int)* pq);

	wmean = mean(x, N);

	for (i = 0; i < N; ++i) {
		inp[i] = x[i] - wmean;
	}

	// Estimate AR(m) coefficients

	ywalg(inp, N, m, phim);

	for (i = 0; i < N; ++i) {
		a[i] = 0.0;
	}

	for (t = m; t < N; ++t) {
		a[t] = inp[t];
		for (i = 0; i < m; ++i) {
			a[t] -= phim[i] * inp[t - i - 1];
		}
	}

	for (i = 0; i < N - m - q; ++i) {
		t = i * pq;
		for (j = 0; j < p; ++j) {
			z[t + j] = inp[m + q + i - j - 1];
		}

		for (k = 0; k < q; ++k) {
			z[t + p + k] = -a[m + q + i - k - 1];
		}
	}

	mtranspose(z, N - m - q, pq, zt);
	mmult(zt, z, A, pq, N - m - q, pq);
	mmult(zt, inp + m + q, b, pq, N - m - q, 1);


	ludecomp(A, pq, ipiv);
	linsolve(A, pq, b, ipiv, temp);

	for (i = 0; i < p; ++i) {
		phi[i] = temp[i];
	}

	for (i = p; i < pq; ++i) {
		theta[i-p] = temp[i];
	}
	sos = 0.0;
	for (t = m+q; t < N; ++t) {
		wmean = inp[t];
		for (i = 0; i < p; ++i) {
			wmean -= (phi[i] * inp[t - i - 1]);
		}

		for (i = 0; i < q; ++i) {
			wmean += (theta[i] * a[t - i - 1]);
		}

		sos += (wmean*wmean);
	}

	*var = sos / (N - m - q);

	free(inp);
	free(phim);
	free(a);
	free(z);
	free(zt);
	free(A);
	free(b);
	free(ipiv);
	free(temp);
}

void hrstep2(double *inp, int N, int m, int p, int q, int pq, double *a, double *sos, double *var) {
	double *z, *zt, *A, *b, *temp;
	int *ipiv;
	int i,j,k,t;
	double wmean;

	z = (double*)malloc(sizeof(double)* (N - m - q) * pq);
	zt = (double*)malloc(sizeof(double)* (N - m - q) * pq);
	A = (double*)malloc(sizeof(double)* pq * pq);
	b = (double*)malloc(sizeof(double)* pq);
	temp = (double*)malloc(sizeof(double)* pq);
	ipiv = (int*)malloc(sizeof(int)* pq);

	for (i = 0; i < N - m - q; ++i) {
		t = i * pq;
		for (j = 0; j < p; ++j) {
			z[t + j] = inp[m + q + i - j - 1];
		}

		for (k = 0; k < q; ++k) {
			z[t + p + k] = -a[m + q + i - k - 1];
		}
	}

	mtranspose(z, N - m - q, pq, zt);
	mmult(zt, z, A, pq, N - m - q, pq);
	mmult(zt, inp + m + q, b, pq, N - m - q, 1);


	ludecomp(A, pq, ipiv);
	linsolve(A, pq, b, ipiv, temp);

	*sos = 0.0;

	for (t = m + q; t < N; ++t) {
		wmean = inp[t];
		for (i = 0; i < p; ++i) {
			wmean -= (temp[i] * inp[t - i - 1]);
		}

		for (i = 0; i < q; ++i) {
			wmean += (temp[i+p] * a[t - i - 1]);
		}

		*sos += (wmean*wmean);
	}

	*var = *sos / (N - m - q);

	free(z);
	free(zt);
	free(A);
	free(b);
	free(ipiv);
	free(temp);
}


void pacf_burg(double* vec, int N, double* par, int M) {
	int i;
	double *temp;
	double var;

	temp = (double*)malloc(sizeof(double)* M);
	for (i = 0; i < M; ++i) {
		burgalg(vec, N, i + 1, temp, &var);
		par[i] = temp[i];
	}
	free(temp);
}

void pacf_yw(double* vec, int N, double* par, int M) {
	int i;
	double *temp;

	temp = (double*)malloc(sizeof(double)* M);
	for (i = 0; i < M; ++i) {
		ywalg(vec, N, i + 1, temp);
		par[i] = temp[i];
	}
	free(temp);
}

void pacf_mle(double* vec, int N, double* par, int M) {
	// Using Box-Jenkins Conditional MLE
	double var, thetac, eps,loglik;
	double *phi, *theta, *varcovar, *res;
	int i, p, mean, k,ret;

	res = (double*)malloc(sizeof(double)* N);
	theta = (double*)malloc(sizeof(double)* 0);

	eps = macheps();
	mean = 1;

	for (i = 0; i < M; ++i) {
		p = i + 1;
		k = p + mean;
		phi = (double*)malloc(sizeof(double)* p);
		varcovar = (double*)malloc(sizeof(double)* k * k);
		ret = nlalsm(vec, N, p, 0, 0, phi, theta, mean, &thetac, &var, eps, varcovar, res);
		//ret = css(vec, N, 6, p, 0, 0, phi, theta, &thetac, &var, res, &loglik, varcovar);
		par[i] = phi[i];
		free(phi);
		free(varcovar);
	}

	free(theta);
	free(res);
}


