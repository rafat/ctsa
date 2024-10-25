// SPDX-License-Identifier: BSD-3-Clause
/*
 * talg.c
 *
 *  Created on: Jul 31, 2013
 *      Author: USER
 */

#include "talg.h"

void detrend_ma(double *sig,int N, int window,double *oup) {
	// Uses moving average filter to detrend the signal

	mafilter(sig,N,window,oup);

	// other options are
	// 1.  mafilter2(sig,N,window,oup)
	// 2.  expfilter(sig,N,window,oup)
	// 3. mafilter_wt(sig,N,weights,window,oup) [See the code in filter.c]
	// For more information, please consult
	//Time Series: Theory and Methods Second Edition
	// by Peter J. Brockwell ,Richard A. Davis
	// Chapter 1
}

int poly(double *A, double *B, double *C, int lA, int lB) {
	int lC,i,j,k;
	double temp;
	lC = lA + lB - 1;
	
	for(i = 0; i < lC;++i) {
			C[i] = 0.;
	}
	
	for(i = 0; i < lC;++i) {
		temp = 0.0;
		for(j = 0; j < lA;++j) {
			for(k = 0; k < lB;++k) {
				if (j + k == i) {
					temp += A[j] * B[k];
				}
			}
		}
		C[i] = temp;
			
	}
	
		
	return lC;
}

int upsample(double *x,int lenx, int M, double *y) {
	int N,i,j,k;

	if (M < 0) {
		return -1;
	}

	if (M == 0) {
		for (i = 0; i < lenx; ++i) {
			y[i] = x[i];
		}
		return lenx;
	}

	N = M * (lenx-1) + 1;
	j = 1;
	k = 0;

	for(i = 0; i < N;++i) {
		j--;
		y[i] = 0.0;
		if (j == 0) {
			y[i] = x[k];
			k++;
			j = M;
		}
	}

	return N;
}

int downsample(double *x,int lenx, int M, double *y) {
	int N,i;

	if (M < 0) {
		return -1;
	}
	if (M == 0) {
		for (i = 0; i < lenx; ++i) {
			y[i] = x[i];
		}
		return lenx;
	}

	N = (lenx-1)/M + 1;

	for (i = 0;i < N;++i) {
		y[i] = x[i*M];
	}

	return N;
}

void deld(int d, double* C) {
	int i,j;
	double *vec,*oup;
	/*
	 * deld returns coefficients of [1,-1]^d
	 * eg.,
	 * d=2 yields [1,-2,1]
	 * d=3 [1,-3,3,1]
	 */ 
	vec = (double*) malloc(sizeof(double) * 2);
	oup = (double*) malloc(sizeof(double) * (d+1));
	vec[0] = 1.;
	vec[1] = -1.;
	oup[0] = 1.;
	// oup has length (d+1)
	
	for(i = 0; i < d; ++i) {
		poly(oup,vec,C,i+1,2);
		for(j = 0; j < i+2; ++j) {
			oup[j] = C[j];
		}
		
	}
	free(vec);
	free(oup);
}

void delds(int D, int s, double *C) {
	int i,j;
	double *vec,*oup;
	/*
	 * delds returns coefficients of [1,....(s-1) zeroes,.. -1]^D
	 * eg., D =1 s = 4 yields
	 * [1,0,0,0,-1]
	 * D=2 s=4 -> [1,0,0,0,-2,0,0,0,1]
	 */ 
	/*
	 * s - seasonal component
	 * s = 4 (quarterly data), 12 (yearly data) etc.
	 * D - Differencing operator. D >= 1 when seasonal component is present
	 * D is usually 1 when seasonal component is present. D = 0 means no
	 * seasonal component
	 */ 
	vec = (double*) malloc(sizeof(double) * (s+1));
	oup = (double*) malloc(sizeof(double) * (D*s + 1));
	for(i = 0; i < s+1; ++i) {
		vec[i] = 0.;
	}
	vec[0] = 1.0;
	vec[s] = -1.0;
	oup[0] = 1.0;
	
	for(i = 0; i < D; ++i) {
		poly(oup,vec,C,i*s+1,s+1);
		for(j = 0; j < s*(i+1)+1; ++j) {
			oup[j] = C[j];
		}
		
	}
	free(vec);
	free(oup);
}

int diff(double *sig, int N, int d, double *oup) {
	int Noup,i,j;
	double *coeff;
	double sum;
	/*
	 * 
	 * diff output = [ X(t) - X(t-1) ] ^ d where X(t) is the
	 * input timeseries of length N
	 * output is of length N - d
	 * 
	 * diff is used to detrend data
	 * 
	 * d=1 typically takes care of linear trend while
	 * d=2 handles quadratic trend
	 */ 
	
	coeff = (double*) malloc(sizeof(double) * (d+1));
	deld(d,coeff);
	Noup = N - d;
	
	for(i = d; i < N;++i) {
		sum = 0.;
		for(j = 1; j < d+1;++j) {
			sum += sig[i-j]*coeff[j];
		}
		oup[i-d] = sum + coeff[0] * sig[i];
		
	}

	free(coeff);
	
	
	return Noup;
}

int diffs(double *sig, int N, int D,int s, double *oup) {
	int Noup,i,j,d;
	double *coeff;
	double sum;
	/*
	 * diffs output = [ X(t) - X(t-s) ] ^ D where X(t) is the
	 * input timeseries of length N
	 * output is of length N - D*s
	 * 
	 * 
	 * s - seasonal component, D = 1 when seasonal component is present and 0
	 * when it is absent. D >= 2 is rarely needed and in most cases you may want
	 *  to re-analyze data if D >= 2 is needed.
	 * s = 4 (quarterly data), 12 (yearly data) etc.
	 * 
	 * D - Differencing operator. D >= 1 when seasonal component is present
	 * D is usually 1 when you are differencing [ X(t) - X(t-s) ]
	 * D = 2 yields [ X(t) - X(t-s) ] * [ X(t) - X(t-s) ]
	 * see Chapter 9 of Box Jenkins for detrending and deseasoning seasonal data
	 * with trends 
	 */ 
	
	d = D*s;
	coeff = (double*) malloc(sizeof(double) * (d+1));
	delds(D,s,coeff);
	Noup = N - d;
	
	for(i = d; i < N;++i) {
		sum = 0.;
		for(j = s; j < d+1;j+=s) {
			sum += sig[i-j]*coeff[j];
		}
		oup[i-d] = sum + coeff[0] * sig[i];
		
	}
	
	free(coeff);
	
	return Noup;
}

void deseason_ma(double *sig,int N,int s,double *oup) {
	double *mt,*w,*seas;
	int k,window,q,odd,jd,count,per,it;
	double temp,wi;
	
	window = s;
	mt = (double*) malloc(sizeof(double) * N);	
	w = (double*) malloc(sizeof(double) * window);	
	seas = (double*) malloc(sizeof(double) * window);
	
	mafilter2(sig,N,window,mt);
	odd = window - ((window/2) * 2);
	
	if (odd) {
		q = (window - 1 ) / 2; 
	} else {
		q = window / 2;
	}
	
	for(k = 0; k < window;++k) {
		jd = k;
		count = 0;
		temp = 0.0;
		while (jd < N - 2*q) {
			temp = temp + sig[q+jd]-mt[q+jd];
			count++;
			jd+=window;
		}
		per = k + q -((k+q)/window)*window;
		w[per] = temp/count;
	}
	
	temp = 0.0;
	for(k = 0; k < window;++k)  {
		temp += w[k];
	}
	wi = temp/window;
	
	for(k = 0; k < window;++k)  {
		seas[k] = w[k] - wi;
	}
	
	for(k = 0; k < N;++k) {
		it = k;
		while (it >= window) {
			it -= window;
		}
		oup[k] = sig[k] - seas[it];
	}
	
	free(w);
	free(seas);
	free(mt);
}

void psiweight(double *phi,double *theta,double *psi,int p,int q,int j) {
	int i,k;
	double temp;
	double *th;
	psi[0] = 1.0;
	th = (double*) malloc(sizeof(double) * (q+1));	
	th[0] = 1.;
	for(i = 0; i < q;++i) {
		th[i+1] = theta[i];
	}
	
	for(i = 1; i < j;++i) {
		psi[i] = 0.0;
		temp = 0.0;
		if(i <= q) {
			psi[i] = th[i];
		}
		for(k = 1; k < p+1;++k) {
			if((i - k) >= 0) {
				temp+=phi[k-1] * psi[i-k];
			} 
			
		}
		psi[i] += temp;
	}
	
	free(th);
}

void piweight(double *phi,double *theta,double *piw,int p,int q,int j) {
	int i,k;
	double temp;
	double *ph;
	piw[0] = 1.0;
	ph = (double*) malloc(sizeof(double) * (p+1));	
	ph[0] = -1.;
	for(i = 0; i < p;++i) {
		ph[i+1] = phi[i];
	}
	
	for(i = 1; i < j;++i) {
		piw[i] = 0.0;
		temp = 0.0;
		if(i <= p) {
			piw[i] = -ph[i];
		}
		for(k = 1; k < q+1;++k) {
			if((i - k) >= 0) {
				temp+=theta[k-1] * piw[i-k];
			} 
			
		}
		piw[i] -= temp;
	}
	
	free(ph);
	
}

void arma_autocovar(double *phi,double *theta,int p,int q,double var,double* acov, int lag) {
	int i,j,t,m;
	double *tcov,*psi,*A,*b,*ph,*th,*thph;
	int *ipiv;
	int p1;
	double temp;
	
	p1 = p+1;
	tcov = (double*) malloc(sizeof(double) * (p1));
	ipiv = (int*) malloc(sizeof(int) * (p1));
	A = (double*) malloc(sizeof(double) * (p1) * (p1));
	b = (double*) malloc(sizeof(double) * (p1));
	ph = (double*) malloc(sizeof(double) * (p1));
	psi = (double*) malloc(sizeof(double) * (q+1));
	th = (double*) malloc(sizeof(double) * (q+1));

	ph[0] = -1.0;
	th[0] = 1.0;
	
	for(i = 0; i < p;++i) {
		ph[i+1] = phi[i];
	}
	
	for(i = 0; i < q;++i) {
		th[i+1] = theta[i];
	}
	
	if (p >= q+1) {
		m = p;
	} else {
		m = q + 1;
	}
	thph = (double*) malloc(sizeof(double) * m);
	// set A
	for(i = 0; i < p1;++i) {
		for(j=0; j < p1;++j) {
			A[i*p1+j] = 0.;
			if (i == j && i != 0) {
				A[i*p1+j] = 1.;
			}

			t = i - j;
			if (t < 0) {
				t = -t;
			} 
			
			A[i*p1+t] -= ph[j];
			
		}
	}
	
	//set b
	
	psiweight(phi,theta,psi,p,q,q+1);
	for(i=0; i < m;++i) {
		
		temp = 0.;
		for(j = 0; j < q+1;++j) {
			if(i+j < q+1) {
				temp += th[i+j] * psi[j];
			}
		}
		thph[i] = temp*var;
		
	}
	
	for(i=0; i < p1;++i) {
		b[i] = thph[i];
	}
	
	ludecomp(A,p1,ipiv);
	linsolve(A,p1,b,ipiv,tcov);
	
	for(i = 0; i < p1;++i) {
		acov[i] = tcov[i];
	}
	
	for(i = p1; i < lag;++i) {
		temp = 0.0;
		for(j = 1; j < p1; ++j) {
			temp+= phi[j-1] * acov[i - j];
		}
		if (i < m) {
			temp += thph[i];
		}
		acov[i] = temp;
	}
	
	free(ph);
	free(th);
	free(thph);
	free(psi);
	free(tcov);
	free(A);
	free(ipiv);
	free(b);
}

int twacf(double *P, int MP, double *Q, int MQ, double *ACF, int MA, double *CVLI, int MXPQ1, double *ALPHA, int MXPQ) {
	int ifault, i, k, kc, j, jpk, kcp1mj, j1, kp1, kp2mj, miim1p, imj, mikp, kp1mj;
	double epsil2, zero, one, half, two, div;

	ifault = 0;
	epsil2 = 1.0e-10;
	zero = 0.0; half = 0.5, one = 1.0, two = 2.0;

	if (MP < 0 || MQ < 0) {
		ifault = 1;
	}
	if (MXPQ != imax(MP, MQ)) {
		ifault = 2;
	}
	if (MXPQ1 != MXPQ + 1) {
		ifault = 3;
	}
	if (MA < MXPQ1) {
		ifault = 4;
	}

	if (ifault > 0) {
		return ifault;
	}

	// Initialization and return if MP = MQ = 0

	ACF[0] = one;
	CVLI[0] = one;

	if (MA == 1) {
		return ifault;
	}

	for (i = 1; i < MA; ++i) {
		ACF[i] = zero;
	}

	if (MXPQ1 == 1) {
		return ifault;
	}

	for (i = 1; i < MXPQ1; ++i) {
		CVLI[i] = zero;
	}
	for (k = 0; k < MXPQ; ++k) {
		ALPHA[k] = 0.0;
	}

	// Computation of the A.C.F. of the moving average part stored in ACF

	if (MQ != 0) {
		for (k = 1; k <= MQ; ++k) {
			CVLI[k] = -Q[k - 1];
			ACF[k] = -Q[k - 1];
			kc = MQ - k;
			if (kc != 0) {
				for (j = 1; j <= kc; ++j) {
					jpk = j + k;
					ACF[k] += (Q[j - 1] * Q[jpk - 1]);
				}
			}//120
			ACF[0] += (Q[k - 1] * Q[k - 1]);
		}

		//Initialization of CVLI = T.W.-S.PHI -- return if MP = 0
	}//180

	if (MP == 0) {
		return ifault;
	}

	for (k = 0; k < MP; ++k) {
		ALPHA[k] = P[k];
		CVLI[k] = P[k];
	}

	// Computation of T.W.-S ALPHA and DELTA
	// DELTA stored in ACF which is gradually overwritten

	for (k = 1; k <= MXPQ; ++k) {
		kc = MXPQ - k;
		if (kc < MP) {
			div = one - ALPHA[kc] * ALPHA[kc];
			if (div <= epsil2) {
				return 5;
			}
			if (kc == 0) {
				break; //break For loop. Go to 290
			}//290
			for (j = 1; j <= kc; ++j) {
				kcp1mj = kc - j;
				ALPHA[j - 1] = (CVLI[j - 1] + ALPHA[kc] * CVLI[kcp1mj]) / div;
			}
		}//240
		if (kc < MQ) {
			j1 = imax(kc + 1 - MP, 1);
			for (j = j1; j <= kc; ++j) {
				kcp1mj = kc - j;
				ACF[j] += ACF[kc + 1] * ALPHA[kcp1mj];
			}
		}//260
		if (kc < MP) {
			for (j = 1; j <= kc; ++j) {
				CVLI[j - 1] = ALPHA[j - 1];
			}
		}
	}//290

	// Computation of T.W.-S NU
	// NU is stored in CVLI copied into ACF

	ACF[0] *= half;
	for (k = 1; k <= MXPQ; ++k) {
		if (k <= MP) {
			kp1 = k + 1;
			div = one - ALPHA[k - 1] * ALPHA[k - 1];
			for (j = 1; j <= kp1; ++j) {
				kp2mj = k + 2 - j;
				CVLI[j - 1] = (ACF[j - 1] + ALPHA[k - 1] * ACF[kp2mj - 1]) / div;
			}
			for (j = 1; j <= kp1; ++j) {
				ACF[j - 1] = CVLI[j - 1];
			}
		}//330
	}//330

	//Computation of ACF

	for (i = 1; i <= MA; ++i) {
		miim1p = imin(i - 1, MP);
		if (miim1p != 0) {
			for (j = 1; j <= miim1p; ++j) {
				imj = i - j; 
				ACF[i - 1] += P[j - 1] * ACF[imj - 1];
			}
		}//430
	}//430

	ACF[0] *= two;

	//Computation of CVLI

	CVLI[0] = one;
	if (MQ > 0) {
		for (k = 1; k <= MQ; ++k) {
			CVLI[k] = -Q[k - 1];
			if (MP != 0) {
				mikp = imin(k, MP);
				for (j = 1; j <= mikp; ++j) {
					kp1mj = k + 1 - j;
					CVLI[k] += P[j - 1] * CVLI[kp1mj - 1];
				}
			}
		}
	}

	return ifault;
}

void artrans(int p, double *old, double *new1) {
	int j, k;
	double a;
	double *temp;

	temp = (double*)malloc(sizeof(double)* p);

	for (j = 0; j < p; ++j) {
		new1[j] = tanh(old[j]);
		temp[j] = new1[j];
	}

	for (j = 1; j < p; ++j) {
		a = new1[j];
		for (k = 0; k < j; ++k) {
			temp[k] -= a * new1[j - k - 1];
		}
		for (k = 0; k < j; ++k) {
			new1[k] = temp[k];
		}
	}

	free(temp);
}

void arinvtrans(int p, double *old, double *new1) {
	int j, k;
	double a;
	double *temp;

	temp = (double*)malloc(sizeof(double)* p);

	for (j = 0; j < p; ++j) {
		temp[j] = new1[j] = old[j];
	}

	for (j = p - 1; j > 0; --j) {
		a = new1[j];
		for (k = 0; k < j; ++k) {
			temp[k] = (new1[k] + a * new1[j - k - 1]) / (1 - a * a);
		}
		for (k = 0; k < j; ++k) {
			new1[k] = temp[k];
		}
	}

	for (j = 0; j < p; ++j) {
		new1[j] = atanh(new1[j]);
	}
	free(temp);
}

void gradtrans(double *raw, int p, int q, int P, int Q,int M, double *A) {
	int i, j, v, N;
	double w1[100], w2[100], w3[100];
	double eps = 1e-03;
	N = p + q + P + Q + M;

	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			if (i == j) {
				A[i + j*N] = (double)i;
			}
		}
	}

	if (p > 0) {
		memcpy(w1, raw, sizeof(double)*p);
		artrans(p, w1, w2);

		for (i = 0; i < p; ++i) {
			w1[i] += eps;
			artrans(p, w1, w3);
			for (j = 0; j < p; ++j) {
				A[i + j*N] = (w3[j] - w2[j]) / eps;
			}
			w1[i] -= eps;
		}
	}

	if (q > 0) {
		v = p;
		memcpy(w1, raw+v, sizeof(double)*q);
		artrans(q, w1, w2);

		for (i = 0; i < q; ++i) {
			w1[i] += eps;
			artrans(q, w1, w3);
			for (j = 0; j < q; ++j) {
				A[i + j*N + v] = (w3[j] - w2[j]) / eps;
			}
			w1[i] -= eps;
		}
	}

	if (P > 0) {
		v = p + q;
		memcpy(w1, raw + v, sizeof(double)*P);
		artrans(P, w1, w2);

		for (i = 0; i < P; ++i) {
			w1[i] += eps;
			artrans(P, w1, w3);
			for (j = 0; j < P; ++j) {
				A[i + (j + v)*N + v] = (w3[j] - w2[j]) / eps;
			}
			w1[i] -= eps;
		}
	}

	if (Q > 0) {
		v = p + q + P;
		memcpy(w1, raw + v, sizeof(double)*Q);
		artrans(Q, w1, w2);

		for (i = 0; i < Q; ++i) {
			w1[i] += eps;
			artrans(Q, w1, w3);
			for (j = 0; j < Q; ++j) {
				A[i + (j + v)*N + v] = (w3[j] - w2[j]) / eps;
			}
			w1[i] -= eps;
		}
	}
}

int archeck(int p, double *ar) {
	int check, N, i;
	double *coeff,*zeror,*zeroi,mod,wmaxv;
	int *which,wmax;

	N = p + 1;
	check = 1;

	coeff = (double*)malloc(sizeof(double)*N);
	zeror = (double*)malloc(sizeof(double)*p);
	zeroi = (double*)malloc(sizeof(double)*p);
	which = (int*)malloc(sizeof(int)*N);

	coeff[0] = 1.0;

	for(i = 0; i < p;++i) {
		coeff[i+1] = -ar[i];
	}

	for(i = 0; i < N;++i) {
		if (coeff[i] != 0) {
			which[i] = i;
		} else {
			which[i] = -1;
		}
	}

	wmax = -1;

	for(i = 0; i < N;++i) {
		if (wmax < which[i]) {
			wmax = which[i];
		}
	}

	if (!wmax){
		free(coeff);
		free(zeror);
		free(zeroi);
		free(which);
		return 1;
	}

	//polyroot(double *coeff, int DEGREE, double *ZEROR, double *ZEROI)

	polyroot(coeff,p,zeror,zeroi);

	for(i = 0; i < p;++i) {
		mod = sqrt(pow(zeror[i],2.0) + pow(zeroi[i],2));
		if (mod < 1.0) {
			free(coeff);
			free(zeror);
			free(zeroi);
			free(which);
			return 0;
		}
	}


	free(coeff);
	free(zeror);
	free(zeroi);
	free(which);
	return check;
}

int invertroot(int q, double *ma) {
	int i, index, retval, fail, rcheck, qn, i1, j;
	double *temp, *zeror, *zeroi, *xr, *xi,*yr,*yi;
	int *ind;
	double mod,tempr,tempi;


	retval = 0;
	index = -1;
	for (i = 0; i < q; ++i) {
		if (ma[i] != 0.0) {
			index = i;
		}
	}

	if (index == -1) {
		return retval;
	}
	
	temp = (double*)malloc(sizeof(double)* (q + 1));
	zeror = (double*)malloc(sizeof(double)* q);
	zeroi = (double*)malloc(sizeof(double)* q);
	ind = (int*)malloc(sizeof(int)* q);

	index++;
	temp[0] = 1.0;
	for (i = 1; i <= index; ++i) {
		temp[i] = ma[i - 1];
	}

	qn = index;

	xr = (double*)malloc(sizeof(double)* (qn + 1));
	xi = (double*)malloc(sizeof(double)* (qn + 1));
	yr = (double*)malloc(sizeof(double)* (qn + 1));
	yi = (double*)malloc(sizeof(double)* (qn + 1));

	fail = polyroot(temp, qn, zeror, zeroi);

	//mdisplay(zeroi, 1, qn);

	if (fail == 1) {
		free(zeror);
		free(zeroi);
		free(temp);
		free(ind);
		free(xr);
		free(xi);
		free(yr);
		free(yi);
		return retval;
	}

	rcheck = 0;

	for (i = 0; i < qn; ++i) {
		mod = zeror[i] * zeror[i] + zeroi[i] * zeroi[i];
		ind[i] = 0;
		if (mod < 1.0) {
			ind[i] = 1;
			rcheck++;
		}
	}

	if (rcheck == 0) {
		free(zeror);
		free(zeroi);
		free(temp);
		free(ind);
		free(xr);
		free(xi);
		free(yr);
		free(yi);
		return retval;
	}
	else {
		retval = 1;
	}

	if (index == 1) {
		ma[0] = 1.0 / ma[0];
		for (i = 1; i < q; ++i) {
			ma[i] = 0.0;
		}
		free(zeror);
		free(zeroi);
		free(temp);
		free(ind);
		free(xr);
		free(xi);
		free(yr);
		free(yi);
		return retval;
	}

	for (i = 0; i < index; ++i) {
		if (ind[i] == 1) {
			mod = zeror[i] * zeror[i] + zeroi[i] * zeroi[i];
			zeror[i] = zeror[i] / mod;
			zeroi[i] = -zeroi[i] / mod;
		}
	}

	//mdisplay(zeroi, 1, qn);

	xr[0] = 1.0; xi[0] = 0.0;
	yr[0] = xr[0];
	yi[0] = xi[0];

	for (i = 0; i < qn; ++i) {
		i1 = i + 1;
		mod = zeror[i] * zeror[i] + zeroi[i] * zeroi[i];
		tempr = zeror[i] / mod;
		tempi = -zeroi[i] / mod;
		xr[i1] = xi[i1] = 0.0;
		for (j = 1; j <= i1; ++j) {
			yr[j] = tempr * xr[j - 1] - tempi * xi[j - 1];
			yi[j] = tempr * xi[j - 1] + tempi * xr[j - 1];
			yr[j] = xr[j] - yr[j];
			yi[j] = xi[j] - yi[j];
		}
		for (j = 1; j <= i1; ++j) {
			xr[j] = yr[j];
			xi[j] = yi[j];
		}
		//mdisplay(xr, 1, qn+1);
	}

	for (i = 0; i < qn; ++i) {
		ma[i] = xr[i+1];
	}

	for (i = qn; i < q; ++i) {
		ma[i] = 0.0;
	}

	free(zeror);
	free(zeroi);
	free(temp);
	free(ind);
	free(xr);
	free(xi);
	free(yr);
	free(yi);
	return retval;
}

void transall(int p,int q, int P, int Q, double *old, double *new1) {
	int N;

	N = 0;

	if (p != 0) {
		artrans(p, old, new1);
		N = N + p;
	}

	if (q != 0) {
		artrans(q, old + N, new1 + N);
		N = N + q;
	}

	if(P != 0) {
		artrans(P, old + N, new1 + N);
		N = N + P;
	}

	if (Q != 0) {
		artrans(Q, old + N, new1 + N);
	}
}

void invtransall(int p, int q, int P, int Q, double *old, double *new1) {
	int N;

	N = 0;

	if (p != 0) {
		arinvtrans(p, old, new1);
		N = N + p;
	}

	if (q != 0) {
		arinvtrans(q, old + N, new1 + N);
		N = N + q;
	}

	if (P != 0) {
		arinvtrans(P, old + N, new1 + N);
		N = N + P;
	}

	if (Q != 0) {
		arinvtrans(Q, old + N, new1 + N);
	}
}

double interpolate_linear(double *xin,double *yin, int N, double z) {
	int i,j,k;
	double *x,*y;
	int *pos;
	double out,ylo,yhi;

	x = (double*)malloc(sizeof(double)*N);
	y = (double*)malloc(sizeof(double)*N);
	pos = (int*)malloc(sizeof(int)*N);

	sort1d_ascending(xin,N,pos);

	for(i = 0; i < N;++i) {
		x[i] = xin[pos[i]];
		y[i] = yin[pos[i]];
	}

	ylo= y[0];
	yhi = y[N-1];

	i = 0;
	j = N - 1;

	if (z < x[0]) {
		free(x);
		free(y);
		free(pos);
		return ylo;
	}

	if (z > x[j]) {
		free(x);
		free(y);
		free(pos);
		return yhi;
	}

	while (i < j-1) {
		k = (i + j)/2;
		if (z < x[k]) {
			j = k;
		} else {
			i = k;
		}
	}

	if (z == x[j]) {
		out = y[j];
		free(x);
		free(y);
		free(pos);
		return out;
	}

	if (z == x[i]) {
		out = y[i];
		free(x);
		free(y);
		free(pos);
		return out;
	}

	out = y[i] + (y[j] - y[i])* ((z - x[i])/(x[j] - x[i]));

	free(x);
	free(y);
	free(pos);

	return out;
}

static double interpolate_linear_sorted_ascending(double *x,double *y, int N, double z) {
	/*
		Input x should be sorted in ascending order with all unique values for x. 
		Input y should have the index ordered on sorted x values.
	*/
	int i,j,k;
	double ylo,yhi;

	i = 0;
	j = N - 1;

	ylo= y[0];
	yhi = y[N-1];

	if (z < x[0]) {
		return ylo;
	}

	if (z > x[j]) {
		return yhi;
	}

	while (i < j-1) {
		k = (i + j)/2;
		if (z < x[k]) {
			j = k;
		} else {
			i = k;
		}
	}

	if (z == x[j]) {
		return y[j];
	}

	if (z == x[i]) {
		return y[i];
	}

	return y[i] + (y[j] - y[i])* ((z - x[i])/(x[j] - x[i]));
}

void linspace(double *x, int N,double xlo,double xhi) {
    int i;
    double intv;

    intv = (xhi - xlo)/(double)(N-1);
    x[0] = xlo;
    x[N-1] = xhi;

    for(i = 1; i < N-1;++i) {
        x[i] = xlo + intv * (double) i;
    }
}

void approx(double *x,double *y, int N,double *xout, double *yout,int Nout) {
    int i;
	/* x and y should be sorted in ascending order. x has all unique values
		xout contains equally spaced Nout values from lowest x[0]
		to highest x[N-1].Outputs y[out] are Nout interpolated values.
	*/

    for(i = 0; i < Nout;++i) {
        if (xout[i] == xout[i]) {
            yout[i] = interpolate_linear_sorted_ascending(x,y,N,xout[i]);
        } else {
            yout[i] = xout[i];
        }
    }

}

void arrayminmax(double *x, int N, double *amin,double *amax) {
	int i;
	*amax = - DBL_MAX;
	*amin = DBL_MAX;

	for(i = 0; i < N;++i) {
		*amax = x[i] > *amax ? x[i] : *amax;
		*amin = x[i] < *amin ? x[i] : *amin;
	}

}

void cumsum(double *x, int N, double *csum) {
	int i;

	double sum = 0.0;

	for (i = 0; i < N;++i) {
		sum += x[i];
		csum[i] = sum;
	}
}

void ppsum(double* u, int n, int l, double* sum)
{
	/* Copyright (C) 1997-2000  Adrian Trapletti 
	efficient computation of the sums involved in the Phillips-Perron tests */
  int i, j;
  double tmp1, tmp2;

  //printf("%d %d %g \n ",n,l,*sum);
  
  tmp1 = 0.0;
  for (i=1; i<= l; i++)
  {
    tmp2 = 0.0;
    for (j=i; j< n; j++)  
    {
      tmp2 += (u[j]*u[j-i]);
    }
    tmp2 *= 1.0-((double)i/((double) l +1.0));
    tmp1 += tmp2;
  }
  tmp1 /= (double) n;
  tmp1 *= 2.0;
  (*sum) += tmp1;
}

/* Common Block Declarations */

struct spans_1_ {
    double spans[3];
};

#define spans_1 (*(struct spans_1_ *) &spans_)

struct consts_1_ {
    double big, sml, eps;
};

#define consts_1 (*(struct consts_1_ *) &consts_)

/* Initialized data */

struct {
    double e_1[3];
    } spans_ = { .05f, .2f, .5f };

struct {
    double e_1[3];
    } consts_ = { 1e20f, 1e-7f, .001f };

static int smooth_(int *n, double *x, double *y, double *w, double *span, int *iper, double *vsmlsq, double *smo, double *acvr)
{
    /* System generated locals */
    int i__1;
    double r__1;
    double d__1;

    /* Local variables */
    double a, h__;
    int i__, j, j0, in, it;
    double xm, ym, wt, sy, fbo, fbw;
    int ibw;
    double var, tmp;
    double xti;
    int out;
    double xto;
    double cvar;
    int jper;

    /* Parameter adjustments */
    --acvr;
    --smo;
    --w;
    --y;
    --x;

    /* Function Body */
    xm = 0.f;
    ym = xm;
    var = ym;
    cvar = var;
    fbw = cvar;
    jper = iabs(*iper);
    ibw = *span * .5 * *n + .5;
    if (ibw < 2) {
	ibw = 2;
    }
    it = (ibw << 1) + 1;
    i__1 = it;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = i__;
	if (jper == 2) {
	    j = i__ - ibw - 1;
	}
	xti = x[j];
	if (j >= 1) {
	    goto L10;
	}
	j = *n + j;
	xti = x[j] - 1.0;
L10:
	wt = w[j];
	fbo = fbw;
	fbw += wt;
	if (fbw > 0.0) {
	    xm = (fbo * xm + wt * xti) / fbw;
	}
	if (fbw > 0.0) {
	    ym = (fbo * ym + wt * y[j]) / fbw;
	}
	tmp = 0.0;
	if (fbo > 0.0) {
	    tmp = fbw * wt * (xti - xm) / fbo;
	}
	var += tmp * (xti - xm);
	cvar += tmp * (y[j] - ym);
/* L20: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	out = j - ibw - 1;
	in = j + ibw;
	if (jper != 2 && (out < 1 || in > *n)) {
	    goto L60;
	}
	if (out >= 1) {
	    goto L30;
	}
	out = *n + out;
	xto = x[out] - 1.0;
	xti = x[in];
	goto L50;
L30:
	if (in <= *n) {
	    goto L40;
	}
	in -= *n;
	xti = x[in] + 1.0;
	xto = x[out];
	goto L50;
L40:
	xto = x[out];
	xti = x[in];
L50:
	wt = w[out];
	fbo = fbw;
	fbw -= wt;
	tmp = 0.0;
	if (fbw > 0.0) {
	    tmp = fbo * wt * (xto - xm) / fbw;
	}
	var -= tmp * (xto - xm);
	cvar -= tmp * (y[out] - ym);
	if (fbw > 0.0) {
	    xm = (fbo * xm - wt * xto) / fbw;
	}
	if (fbw > 0.0) {
	    ym = (fbo * ym - wt * y[out]) / fbw;
	}
	wt = w[in];
	fbo = fbw;
	fbw += wt;
	if (fbw > 0.0) {
	    xm = (fbo * xm + wt * xti) / fbw;
	}
	if (fbw > 0.0) {
	    ym = (fbo * ym + wt * y[in]) / fbw;
	}
	tmp = 0.0;
	if (fbo > 0.0) {
	    tmp = fbw * wt * (xti - xm) / fbo;
	}
	var += tmp * (xti - xm);
	cvar += tmp * (y[in] - ym);
L60:
	a = 0.0;
	if (var > *vsmlsq) {
	    a = cvar / var;
	}
	smo[j] = a * (x[j] - xm) + ym;
	if (*iper <= 0) {
	    goto L80;
	}
	h__ = 0.0;
	if (fbw > 0.0) {
	    h__ = 1.0 / fbw;
	}
	if (var > *vsmlsq) {
/* Computing 2nd power */
	    d__1 = x[j] - xm;
	    h__ += d__1 * d__1 / var;
	}
	acvr[j] = 0.0;
	a = 1.0 - w[j] * h__;
	if (a <= 0.0) {
	    goto L70;
	}
	acvr[j] = (r__1 = y[j] - smo[j], fabs(r__1)) / a;
	goto L80;
L70:
	if (j <= 1) {
	    goto L80;
	}
	acvr[j] = acvr[j - 1];
L80:
	;
    }
    j = 1;
L90:
    j0 = j;
    sy = smo[j] * w[j];
    fbw = w[j];
    if (j >= *n) {
	goto L110;
    }
L100:
    if (x[j + 1] > x[j]) {
	goto L110;
    }
    ++j;
    sy += w[j] * smo[j];
    fbw += w[j];
    if (j < *n) {
	goto L100;
    }
L110:
    if (j <= j0) {
	goto L130;
    }
    a = 0.0;
    if (fbw > 0.0) {
	a = sy / fbw;
    }
    i__1 = j;
    for (i__ = j0; i__ <= i__1; ++i__) {
	smo[i__] = a;
/* L120: */
    }
L130:
    ++j;
    if (j <= *n) {
	goto L90;
    }
    return 0;
} /* smooth_ */


static int supsmu_(int *n, double *x, double *y, double *w, int *iper, double *span, double *alpha, double *smo, double *sc)
{
    /* System generated locals */
    int sc_dim1, sc_offset, i__1;
    double r__1, r__2;
    double d__1, d__2;

    /* Local variables */
    double a, f, h__;
    int i__, j;
    double sw, sy;
    int jper;
    double scale, resmin;
    double vsmlsq;


/* ------------------------------------------------------------------ */

/* super-smoother. */

/* Friedman J.H. (1984). A variable span smoother. Department of Statistics, */
/*    Stanford University, Technical Report LCS5. */

/* version 10/10/84. */

/* coded  and copyright (c) 1984 by: */

/*                        Jerome H. Friedman */
/*                     Department of Statistics */
/*                               and */
/*                Stanford Linear Accelerator Center */
/*                        Stanford University */

/* all rights reserved. */


/* input: */
/*    n : number of observations (x,y - pairs). */
/*    x(n) : ordered abscissa values. */
/*    y(n) : corresponding ordinate (response) values. */
/*    w(n) : weight for each (x,y) observation. */
/*    iper : periodic variable flag. */
/*       iper=1 => x is ordered interval variable. */
/*       iper=2 => x is a periodic variable with values */
/*                 in the range (0.0,1.0) and period 1.0. */
/*    span : smoother span (fraction of observations in window). */
/*           span=0.0 => automatic (variable) span selection. */
/*    alpha : controles high frequency (small span) penality */
/*            used with automatic span selection (bass tone control). */
/*            (alpha.le.0.0 or alpha.gt.10.0 => no effect.) */
/* output: */
/*   smo(n) : smoothed ordinate (response) values. */
/* scratch: */
/*   sc(n,7) : internal working storage. */

/* note: */
/*    for small samples (n < 40) or if there are substantial serial */
/*    correlations between obserations close in x - value, then */
/*    a prespecified fixed span smoother (span > 0) should be */
/*    used. reasonable span values are 0.2 to 0.4. */

/* ------------------------------------------------------------------ */

    /* Parameter adjustments */
    sc_dim1 = *n;
    sc_offset = 1 + sc_dim1;
    sc -= sc_offset;
    --smo;
    --w;
    --y;
    --x;

    /* Function Body */
    if (x[*n] > x[1]) {
	goto L30;
    }
    sy = 0.0;
    sw = sy;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sy += w[j] * y[j];
	sw += w[j];
/* L10: */
    }
    a = 0.0;
    if (sw > 0.0) {
	a = sy / sw;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	smo[j] = a;
/* L20: */
    }
    return 0;
L30:
    i__ = *n / 4;
    j = i__ * 3;
    scale = x[j] - x[i__];
L40:
    if (scale > 0.f) {
	goto L50;
    }
    if (j < *n) {
	++j;
    }
    if (i__ > 1) {
	--i__;
    }
    scale = x[j] - x[i__];
    goto L40;
L50:
/* Computing 2nd power */
    r__1 = consts_1.eps * scale;
    vsmlsq = r__1 * r__1;
    jper = *iper;
    if (*iper == 2 && (x[1] < 0.0 || x[*n] > 1.0)) {
	jper = 1;
    }
    if (jper < 1 || jper > 2) {
	jper = 1;
    }
    if (*span <= 0.0) {
	goto L60;
    }
    smooth_(n, &x[1], &y[1], &w[1], span, &jper, &vsmlsq, &smo[1], &sc[
	    sc_offset]);
    return 0;
L60:
    for (i__ = 1; i__ <= 3; ++i__) {
	smooth_(n, &x[1], &y[1], &w[1], &spans_1.spans[i__ - 1], &jper, &
		vsmlsq, &sc[((i__ << 1) - 1) * sc_dim1 + 1], &sc[sc_dim1 * 7 
		+ 1]);
	i__1 = -jper;
	smooth_(n, &x[1], &sc[sc_dim1 * 7 + 1], &w[1], &spans_1.spans[1], &
		i__1, &vsmlsq, &sc[(i__ << 1) * sc_dim1 + 1], &h__);
/* L70: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	resmin = consts_1.big;
	for (i__ = 1; i__ <= 3; ++i__) {
	    if (sc[j + (i__ << 1) * sc_dim1] >= resmin) {
		goto L80;
	    }
	    resmin = sc[j + (i__ << 1) * sc_dim1];
	    sc[j + sc_dim1 * 7] = spans_1.spans[i__ - 1];
L80:
	    ;
	}
	if (*alpha > 0.f && *alpha <= 10.f && resmin < sc[j + sc_dim1 * 6] && 
		resmin > 0.f) {
/* Computing MAX */
	    r__1 = consts_1.sml, r__2 = resmin / sc[j + sc_dim1 * 6];
	    d__1 = (double) pmax(r__1,r__2);
	    d__2 = (double) (10.0 - *alpha);
	    sc[j + sc_dim1 * 7] += (spans_1.spans[2] - sc[j + sc_dim1 * 7]) * 
		    pow(d__1, d__2);
	}
/* L90: */
    }
    i__1 = -jper;
    smooth_(n, &x[1], &sc[sc_dim1 * 7 + 1], &w[1], &spans_1.spans[1], &i__1, &
	    vsmlsq, &sc[(sc_dim1 << 1) + 1], &h__);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (sc[j + (sc_dim1 << 1)] <= spans_1.spans[0]) {
	    sc[j + (sc_dim1 << 1)] = spans_1.spans[0];
	}
	if (sc[j + (sc_dim1 << 1)] >= spans_1.spans[2]) {
	    sc[j + (sc_dim1 << 1)] = spans_1.spans[2];
	}
	f = sc[j + (sc_dim1 << 1)] - spans_1.spans[1];
	if (f >= 0.f) {
	    goto L100;
	}
	f = -f / (spans_1.spans[1] - spans_1.spans[0]);
	sc[j + (sc_dim1 << 2)] = (1.f - f) * sc[j + sc_dim1 * 3] + f * sc[j + 
		sc_dim1];
	goto L110;
L100:
	f /= spans_1.spans[2] - spans_1.spans[1];
	sc[j + (sc_dim1 << 2)] = (1.f - f) * sc[j + sc_dim1 * 3] + f * sc[j + 
		sc_dim1 * 5];
L110:
	;
    }
    i__1 = -jper;
    smooth_(n, &x[1], &sc[(sc_dim1 << 2) + 1], &w[1], spans_1.spans, &i__1, &
	    vsmlsq, &smo[1], &h__);
    return 0;
} /* supsmu_ */


/* --------------------------------------------------------------- */

/* this sets the compile time (default) values for various */
/* internal parameters : */

/* spans : span values for the three running linear smoothers. */
/* spans(1) : tweeter span. */
/* spans(2) : midrange span. */
/* spans(3) : woofer span. */
/* (these span values should be changed only with care.) */
/* big : a large representable floating point number. */
/* sml : a small number. should be set so that (sml)**(10.0) does */
/*       not cause floating point underflow. */
/* eps : used to numerically stabilize slope calculations for */
/*       running linear fits. */

/* these parameter values can be changed by declaring the */
/* relevant labeled common in the main program and resetting */
/* them with executable statements. */

/* ----------------------------------------------------------------- */

void supsmu(double *x, int N, double *y,double *w, int periodic,double span, double alpha,double *oup) {
	/* input: */
	/*    N : number of observations (x,y - pairs). */
	/*    x(N) : ordered abscissa values. */
	/*    y(N) : corresponding ordinate (response) values. */
	/*    w(N) : weight for each (x,y) observation. set to NULL if all weights are equal. */
	/*    periodic : periodic variable flag. */
	/*       iper=1 => x is ordered interval variable. */
	/*       iper=2 => x is a periodic variable with values */
	/*                 in the range (0.0,1.0) and period 1.0. */
	/*    span : smoother span (fraction of observations in window). */
	/*           span=0.0 => automatic (variable) span selection. */
	/*    alpha : controles high frequency (small span) penality */
	/*            used with automatic span selection (bass tone control). */
	/*            (alpha.le.0.0 or alpha.gt.10.0 => no effect.). Set to -1.0 */
	/* output: */
	/*   oup(N) : smoothed ordinate (response) values. */
	double *sc,*weights;
	int iper,i;

	sc = (double*)calloc(N*7,sizeof(double));

	if (periodic == 1 || periodic == 2) {
		iper = periodic;
	} else {
		printf("periodic only takes two values - 1 or 2 \n");
		printf(" 1 : x is ordered interval variable \n");
		printf("2 : x is a periodic variable with values in the range (0.0,1.0) and period 1.0.\n");
		exit(-1);
	}

	if (span < 0.0 || span > 1.0) {
		printf("The fractional data span to use (0 <= span < 1).If specified, then use the given fixed span.\n");
        printf("If set to zero (default) then use a variable span.\n");
		exit(-1);
	}

	if (w == NULL) {
		weights = (double*)calloc(N,sizeof(double));
		for(i = 0; i < N;++i) {
			weights[i] = 1.0;
		}
		supsmu_(&N, x, y, weights, &iper, &span, &alpha, oup, sc);
		free(weights);
	} else {
		supsmu_(&N, x, y, w, &iper, &span, &alpha, oup, sc);
	}

	

	free(sc);
}
