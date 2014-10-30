/*
 * boxjenkins.c
 *
 *  Created on: Jun 29, 2014
 *      Author: Rafat Hussain
 */

#include "boxjenkins.h"

void USPE(double *inp,int N,int p, int q, double *phi,double *theta,double *thetac,double *var) {
	int K,i,j,l,t,temp,it1;
	double wmean,tempd,epsilon;
	double *b,*c,*A,*c2,*x,*f,*tau,*T,*h;
	int *ipiv,*ipiv2;

	K = p + q + 1;

	c = (double *) malloc(sizeof(double) * K);
	A = (double *) malloc(sizeof(double) * p * p);
	b = (double *) malloc(sizeof(double) * p);
	ipiv = (int *) malloc(sizeof(int) * p);
	c2 = (double *) malloc(sizeof(double) * (q+1));
	x = (double *) malloc(sizeof(double) * (p+1));
	f = (double *) malloc(sizeof(double) * (q+1));
	tau = (double *) malloc(sizeof(double) * (q+1));
	h = (double *) malloc(sizeof(double) * (q+1));
	T = (double *) malloc(sizeof(double) * (q+1) * (q+1));
	ipiv2 = (int *) malloc(sizeof(int) * (q+1));


	autocovar(inp,N,c,K);
	wmean = mean(inp,N);
	x[0] = -1.0;

	for(i = 0; i < p;++i) {
		phi[i] = 0.0;
	}

	for(i = 0; i < q;++i) {
		theta[i] = 0.0;
	}

	//Step 1
	if (p > 0) {

		for(i = 1; i <= p;++i) {
			t = (i-1) * p;
			for(j = 1; j <= p;++j) {
				temp = q+i-j;
				if (temp < 0) {
					temp = -temp;
				}
				A[t+j-1] = c[temp];
			}
		}

		for(i = 0; i < p;++i) {
			b[i] = c[q+i+1];
		}

		ludecomp(A,p,ipiv);
		linsolve(A,p,b,ipiv,phi);

	}


	for(i = 0; i < p;++i) {
		x[i+1] = phi[i];
	}

	// Step 2
	if (p == 0) {
		for(i = 0; i < q + 1; ++i) {
			c2[i] = c[i];
		}
	} else {
		for(j = 0; j < q + 1;++j) {
			tempd = 0.0;
			for(i = 0; i < p + 1;++i) {
				for(l = 0; l < p + 1; ++l) {
					temp = j + i - l;
					if (temp < 0) {
						temp = -temp;
					}
					tempd += x[i]*x[l]*c[temp];
				}
			}
			c2[j] = tempd;
		}

	}

	//Step 3 Newton Raphson Method

	epsilon = 0.0001;

	for(i = 0; i <= q;++i) {
		tau[i] = 0.0;
		f[i] = 0.0;
	}
	tau[0] = sqrt(c2[0]);

	for(j = 0; j <= q;++j) {
		for(i = 0; i < q-j+1;++i) {
			f[j] += tau[i] * tau[i+j];
		}
		f[j] -= c2[j];
	}

	it1 = 0;

	while (array_max_abs(f,q+1) > epsilon && it1 < 200) {
		it1++;
		for(i = 0; i < q+1;++i) {
			t = i * (q+1);
			for(j = 0; j < q+1;++j) {
				tempd = 0.0;
				if (i + j <= q) {
					tempd=tempd+tau[i+j];
				}
				if (j - i >= 0 ) {
					tempd=tempd+tau[j-i];
				}
				T[t+j] = tempd;
			}
		}
		//mdisplay(T,q+1,q+1);
		ludecomp(T,q+1,ipiv2);
		linsolve(T,q+1,f,ipiv2,h);

		for(i = 0; i < q+1;++i) {
			tau[i] = tau[i] - h[i];
			f[i] = 0.0;
		}

		for(j = 0; j < q+1;++j) {
			for(i = 0; i < q-j+1;++i) {
				f[j]=f[j]+tau[i]*tau[i+j];
			}
			f[j] = f[j] - c2[j];
		}
	}

	*thetac = 0.0;

	if (p > 0) {
		tempd = 0.0;
		for(i = 0; i < p;++i) {
			tempd += x[i+1];
		}
		*thetac = wmean * (1.0 - tempd);
	} else {
		*thetac = wmean;
	}

	for(j = 1; j < q+1;++j) {
		theta[j-1]=-1.0*tau[j]/tau[0];
	}

	*var = tau[0] * tau[0];

	if (q == 0) {
		tempd = 0.0;
		for(i = 0; i < p;++i) {
			tempd = tempd + c[i+1] * x[i+1];
		}
		*var = c[0] - tempd;
	}

	free(c);
	free(A);
	free(b);
	free(ipiv);
	free(c2);
	free(x);
	free(f);
	free(tau);
	free(T);
	free(ipiv2);
	free(h);
}

void USPE_seasonal(double *inp,int N,int s,int p, int q, double *phi,double *theta) {
	int K,i,j,l,t,temp,it1;
	double tempd,epsilon;
	double *b,*c,*A,*c2,*x,*f,*tau,*T,*h,*cv;
	int *ipiv,*ipiv2;

	K = (p + q + 1) * s;

	cv = (double *) malloc(sizeof(double) * K);
	c = (double *) malloc(sizeof(double) * (p+q+1));
	A = (double *) malloc(sizeof(double) * p * p);
	b = (double *) malloc(sizeof(double) * p);
	ipiv = (int *) malloc(sizeof(int) * p);
	c2 = (double *) malloc(sizeof(double) * (q+1));
	x = (double *) malloc(sizeof(double) * (p+1));
	f = (double *) malloc(sizeof(double) * (q+1));
	tau = (double *) malloc(sizeof(double) * (q+1));
	h = (double *) malloc(sizeof(double) * (q+1));
	T = (double *) malloc(sizeof(double) * (q+1) * (q+1));
	ipiv2 = (int *) malloc(sizeof(int) * (q+1));


	autocovar(inp,N,cv,K);
	for(i = 0; i < p+q+1;++i) {
		c[i] = cv[i*s];
	}

	x[0] = -1.0;

	for(i = 0; i < p;++i) {
		phi[i] = 0.0;
	}

	for(i = 0; i < q;++i) {
		theta[i] = 0.0;
	}

	//Step 1
	if (p > 0) {

		for(i = 1; i <= p;++i) {
			t = (i-1) * p;
			for(j = 1; j <= p;++j) {
				temp = q+i-j;
				if (temp < 0) {
					temp = -temp;
				}
				A[t+j-1] = c[temp];
			}
		}

		for(i = 0; i < p;++i) {
			b[i] = c[q+i+1];
		}

		ludecomp(A,p,ipiv);
		linsolve(A,p,b,ipiv,phi);

	}


	for(i = 0; i < p;++i) {
		x[i+1] = phi[i];
	}

	// Step 2
	if (p == 0) {
		for(i = 0; i < q + 1; ++i) {
			c2[i] = c[i];
		}
	} else {
		for(j = 0; j < q + 1;++j) {
			tempd = 0.0;
			for(i = 0; i < p + 1;++i) {
				for(l = 0; l < p + 1; ++l) {
					temp = j + i - l;
					if (temp < 0) {
						temp = -temp;
					}
					tempd += x[i]*x[l]*c[temp];
				}
			}
			c2[j] = tempd;
		}

	}

	//Step 3 Newton Raphson Method

	epsilon = 0.0001;

	for(i = 0; i <= q;++i) {
		tau[i] = 0.0;
		f[i] = 0.0;
	}
	tau[0] = sqrt(c2[0]);

	for(j = 0; j <= q;++j) {
		for(i = 0; i < q-j+1;++i) {
			f[j] += tau[i] * tau[i+j];
		}
		f[j] -= c2[j];
	}

	it1 = 0;

	while (array_max_abs(f,q+1) > epsilon && it1 < 200) {
		it1++;
		for(i = 0; i < q+1;++i) {
			t = i * (q+1);
			for(j = 0; j < q+1;++j) {
				tempd = 0.0;
				if (i + j <= q) {
					tempd=tempd+tau[i+j];
				}
				if (j - i >= 0 ) {
					tempd=tempd+tau[j-i];
				}
				T[t+j] = tempd;
			}
		}
		//mdisplay(T,q+1,q+1);
		ludecomp(T,q+1,ipiv2);
		linsolve(T,q+1,f,ipiv2,h);

		for(i = 0; i < q+1;++i) {
			tau[i] = tau[i] - h[i];
			f[i] = 0.0;
		}

		for(j = 0; j < q+1;++j) {
			for(i = 0; i < q-j+1;++i) {
				f[j]=f[j]+tau[i]*tau[i+j];
			}
			f[j] = f[j] - c2[j];
		}
	}

	for(j = 1; j < q+1;++j) {
		theta[j-1]=-1.0*tau[j]/tau[0];
	}


	free(c);
	free(cv);
	free(A);
	free(b);
	free(ipiv);
	free(c2);
	free(x);
	free(f);
	free(tau);
	free(T);
	free(ipiv2);
	free(h);

}

static double sosqm(double *vec, int N) {
	double sos;
	int i;
	sos = 0.0;

	for (i = 0; i < N;++i) {
		sos += (vec[i] * vec[i]);
	}

	return sos;
}

void avaluem(double *w,int N,int p, int q, double *phi,double *theta,int tval,double *a) {
	int k,i,t,t2,newl,iter;
	double temp;
	double *e,*w2,*alpha;


	e = (double*) malloc(sizeof(double) * N);
	w2 = (double*) malloc(sizeof(double) * tval);

	if (p > q) {
		k = p;
	} else {
		k = q;
	}

	alpha = (double*) malloc(sizeof(double) * (tval+N+k));


	for(i = 0; i < N;++i) {
		e[i] = 0.0;
	}

	for(i = 0; i < N+tval+k;++i) {
		alpha[i] = 0.0;
	}

	for(i = 0; i < N+tval-1;++i) {
		a[i] = 0.0;
	}

	for(i = 0; i < tval;++i) {
		w2[i] = 0.0;
	}

	for(i = N-k-1; i >= 0; --i) {
		temp = 0.0;

		for (t = 0; t < p;++t) {
			temp = temp - phi[t] * w[i+t+1];
		}

		for ( t2 = 0; t2 < q;++t2) {
			temp = temp + theta[t2] * e[t2+i+1];
		}

		e[i] = w[i] + temp;
	}
/*
	for(i = 0; i < N;++i) {
		printf("%g \n",e[i]);
	}
*/

	for(i = 0; i < tval;++i) {
		temp = 0.0;
		for(t = 0; t < p;++t) {
			if (t-i < 0) {
				temp = temp - phi[t]*w2[i-t-1];
			} else {
				temp = temp - phi[t]*w[t-i];
			}
			//printf("temp %g \n",temp);

		}

		for(t2 = 0; t2 < q;++t2) {
			if (t2-i >= 0) {
				temp=temp+theta[t2]*e[t2-i];
			}
		}
		w2[i] = -temp;
	}
/*
	for(i = 0; i < tval;++i) {
		printf("%g \n",w2[i]);
	}
*/


	newl = N + tval;
	

	for (iter = k; iter < newl+k;iter++) {
		temp = 0.0;
		for (t = 0; t < p;++t) {
			if (iter-t-1 >= tval+k ) {
				temp = temp - phi[t] * w[iter-t-1-tval-k];
			} else if (iter-t-1 < tval+k && iter-t-1 >= k) {
				temp = temp - phi[t] * w2[k-iter+tval+t];
			}
		}
		
		for (t2 = 0; t2 < q; ++t2) {
			if (iter-t2-1 >= k ) {
				temp = temp + theta[t2] * alpha[iter-t2-1];
			}
		}

		//printf("%g \n",temp);

		if (iter >= tval+k) {
			alpha[iter] = w[iter-tval-k] + temp;
		} else if (iter >= k){
			alpha[iter] = w2[k-iter+tval-1] + temp;
		}
	}

	for (iter = k+1; iter < newl+k;iter++) {
		a[iter-k-1] = alpha[iter];
	}

	free(e);
	free(w2);
	free(alpha);
}

void avalues(double *w,int N,int p, int q, double *phi,double *theta,int s,int ps,int qs,
		double *phis, double *thetas,int tval,double *a) {

	int k,i,t,t2,newl,iter;
	double temp;
	double *e,*w2,*alpha,*alpha2;


	e = (double*) malloc(sizeof(double) * N);
	w2 = (double*) malloc(sizeof(double) * tval);

	if (p > q) {
		k = p;
	} else {
		k = q;
	}

	if ( k < ps) {
		k = ps;
	}

	if (k < qs) {
		k = qs;
	}

	alpha = (double*) malloc(sizeof(double) * (tval+N+k));
	alpha2 = (double*) malloc(sizeof(double) * (tval+N+k));

	for(i = 0; i < N;++i) {
		e[i] = 0.0;
	}

	for(i = 0; i < N+tval+k;++i) {
		alpha[i] = 0.0;
	}

	for(i = 0; i < N+tval-1;++i) {
		a[i] = 0.0;
	}

	for(i = 0; i < tval;++i) {
		w2[i] = 0.0;
	}

	for(i = N-k-1; i >= 0; --i) {
		temp = 0.0;

		for (t = 0; t < p;++t) {
			temp = temp - phi[t] * w[i+t+1];
		}

		for ( t2 = 0; t2 < q;++t2) {
			temp = temp + theta[t2] * e[t2+i+1];
		}

		e[i] = w[i] + temp;
	}

	for(i = 0; i < tval;++i) {
		temp = 0.0;
		for(t = 0; t < p;++t) {
			if (t-i < 0) {
				temp = temp - phi[t]*w2[i-t-1];
			} else {
				temp = temp - phi[t]*w[t-i];
			}
			//printf("temp %g \n",temp);

		}

		for(t2 = 0; t2 < q;++t2) {
			if (t2-i >= 0) {
				temp=temp+theta[t2]*e[t2-i];
			}
		}
		w2[i] = -temp;
	}

	newl = N + tval;


	for (iter = k; iter < newl+k;iter++) {
		temp = 0.0;
		for (t = 0; t < p;++t) {
			if (iter-t-1 >= tval+k ) {
				temp = temp - phi[t] * w[iter-t-1-tval-k];
			} else if (iter-t-1 < tval+k && iter-t-1 >= k) {
				temp = temp - phi[t] * w2[k-iter+tval+t];
			}
		}

		for (t2 = 0; t2 < q; ++t2) {
			if (iter-t2-1 >= k ) {
				temp = temp + theta[t2] * alpha[iter-t2-1];
			}
		}

		//printf("%g \n",temp);

		if (iter >= tval+k) {
			alpha[iter] = w[iter-tval-k] + temp;
		} else if (iter >= k){
			alpha[iter] = w2[k-iter+tval-1] + temp;
		}
	}
/*
	for (iter = k+1; iter < newl+k;iter++) {
		a[iter-k-1] = alpha[iter];
	}
*/
	for(i = 0; i < newl+k;++i) {
		alpha2[iter] = 0.0;
	}

	for (iter = k; iter < newl+k;iter++) {
		temp = 0.0;

		for (t = 0; t < p;++t) {
			temp = temp - phis[t] * alpha[iter - (t+1) * s];
		}

		for (t2 = 0; t2 < q; ++t2) {
			temp = temp + thetas[t2] * alpha2[iter-(t2+1)*s];
		}

		alpha2[iter] = alpha[iter] + temp;

	}

	for (iter = k+1; iter < newl+k;iter++) {
		a[iter-k-1] = alpha2[iter];
	}


	free(e);
	free(w2);
	free(alpha);
	free(alpha2);

}

int nlalsm(double *vec,int N,int p,int delta, int q, double *phi,double *theta,
		int M,double *thetac,double *var,double eps,double *varcovar,double *residuals) {
	int retval,i,j,l,k,t,t2,it,Norig,qd,itmax;
	double lam,vmul,e1,di,eps2,stepsize,temp,tau,hmax,suminc,sumb;
	double *inp,*b,*at,*at0,*x,*phi0,*theta0,*orig,*A,*g,*D,*h,*binc,*Atemp;
	double *phiinc,*thetainc,*atinc,*resat,*origval;
	int *ipiv;


	eps2 = sqrt(eps);
	lam = 0.01;
	vmul = 2.0;
	e1 = 1.0e-08;
	it = 0;
	Norig = N;
	k = p + q;
	di = eps2;//0.001;//eps2;//0.001
	qd = q + 20 * p;

	if (qd < 15) {
		qd = 15;
	}
	hmax = 1;
	itmax = 200;

	inp = (double*) malloc(sizeof(double) * (N - delta));
	origval = (double*) malloc(sizeof(double) * (N - delta));
	orig = (double*) malloc(sizeof(double) * (N - delta));
	phi0 = (double*) malloc(sizeof(double) * p);
	theta0 = (double*) malloc(sizeof(double) * q);
	phiinc = (double*) malloc(sizeof(double) * p);
	thetainc = (double*) malloc(sizeof(double) * q);
	
	if (delta > 0) {
		N = diff(vec,Norig,delta,inp);
		for (i = 0; i < N;++i) {
			origval[i] = inp[i];
		}
	} else {
		for(i = 0; i < N;++i) {
			inp[i] = vec[i];
			origval[i] = inp[i];
		}
	}
	

	USPE(inp,N,p,q,phi,theta,thetac,var); // Initial Parameters

	if (M == 1) {
		for(i = 0; i < N;++i) {
			inp[i] -= *thetac;
		}
		M = 1;
		k++;
	} else {
		M = 0;
	}
	
	b = (double*) malloc(sizeof(double) * k);
	at = (double*) malloc(sizeof(double) * (N+qd-1));
	at0 = (double*) malloc(sizeof(double) * (N+qd-1));
	atinc = (double*) malloc(sizeof(double) * (N+qd-1));
	resat = (double*) malloc(sizeof(double) * (N+qd-1));
	x = (double*) malloc(sizeof(double) * k * (N+qd-1));
	A = (double*) malloc(sizeof(double) * k * k);
	Atemp = (double*) malloc(sizeof(double) * k * k);
	g = (double*) malloc(sizeof(double) * k);
	D = (double*) malloc(sizeof(double) * k);
	h = (double*) malloc(sizeof(double) * k);
	binc = (double*) malloc(sizeof(double) * k);
	ipiv = (int*) malloc(sizeof(int) * k);
	
	if (M == 1) {
		b[0] = *thetac;
	}
	
	for(i = 0; i < p;++i) {
		b[M+i] = phi[i];
	}
	
	for(i = 0; i < q;++i) {
		b[M+p+i] = theta[i];
	}
	
	//mdisplay(b,1,k);

    avaluem(inp,N,p,q,phi,theta,qd,at);
	
	for(i = 0; i < N+qd-1;++i) {
		resat[i] = at[i];
	}
	//mdisplay(at,1,N+qd-1);
	for(i = 0; i < k;++i) {
		stepsize = b[i]*di;//di;//b[i] * di;
		temp = b[i];
		b[i] += stepsize;
		stepsize = b[i] - temp;
		if (M == 1) {
			for(j = 0;j < p;++j) {
				phi0[j] = b[1+j];
			}
			for(j = 0; j < q;++j) {
				theta0[j] = b[p+1+j];
			}
			for(j = 0; j < N;++j) {
				orig[j] = origval[j] - b[0];
			}
		} else {
			for(j = 0;j < p;++j) {
				phi0[j] = b[j];
			}
			for(j = 0; j < q;++j) {
				theta0[j] = b[p+j];
			}
			for(j = 0; j < N;++j) {
				orig[j] = inp[j];
			}
		}
	    avaluem(orig,N,p,q,phi0,theta0,qd,at0);	
	    //mdisplay(inp,1,N);

		t = i * (N + qd - 1);
		
		for(j = 0; j < N + qd - 1;++j) {
			x[t+j] = (at[j] - at0[j]) / stepsize;
		}
		
		b[i] -= stepsize;
	}

	for(i = 0; i < k;++i) {
		t = i * (N + qd - 1);
		for(j = 0; j < k;++j) {
			t2 = j * (N + qd - 1);
			A[i*k+j] = 0.0;
			for(l = 0; l < N + qd -1;++l) {
				A[i*k+j] += x[t+l] * x[t2+l];
			}

		}
	}

	for(i = 0; i < k;++i) {
		t = i * (N + qd - 1);
		g[i] = 0.0;
		for(l = 0; l < N + qd - 1;++l) {
			g[i] += x[t+l] * at[l];
		}
		D[i] = sqrt(A[i*k+i]);
	}

	tau = array_max(D,k);
	//printf("TAU %g \n",tau);
	//mdisplay(x,k,N+qd-1);

	while (hmax >= e1 && it < itmax) {
		it++;
		for(i = 0; i < k;++i) {
			for(j = 0; j < k;++j) {
				if (it == 1) {
					A[i*k+j] = A[i*k+j] / (D[i] * D[j]);
				}
				if (i == j) {
					A[i*k+j] = 1.0 + lam;
				}
			}
		}

		if (it == 1) {
			for(i = 0; i < k; ++i) {
				g[i] = g[i] / D[i];
			}
		}

		for(i = 0; i < k*k;++i) {
			Atemp[i] = A[i];
		}
		ludecomp(A,k,ipiv);
		linsolve(A,k,g,ipiv,h);
		for(i = 0; i < k*k;++i) {
			A[i] = Atemp[i];
		}

		for(j = 0; j < k;++j) {
			h[j] = h[j] / D[j];
			binc[j] = b[j] + h[j];
		}

		if (M == 1) {
			for(j = 0;j < p;++j) {
				phiinc[j] = binc[1+j];
			}
			for(j = 0; j < q;++j) {
				thetainc[j] = binc[p+1+j];
			}

			for(j = 0; j < N;++j) {
				orig[j] = origval[j] - binc[0];
			}
		} else {
			for(j = 0;j < p;++j) {
				phiinc[j] = binc[j];
			}
			for(j = 0; j < q;++j) {
				thetainc[j] = binc[p+j];
			}
			for(j = 0; j < N;++j) {
				orig[j] = inp[j];
			}
		}

	    avaluem(orig,N,p,q,phiinc,thetainc,qd,atinc);
	    for(i = 0; i < N+qd-1;++i) {
	    	resat[i] = atinc[i];
	    }

	    suminc = sosqm(atinc,N+qd-1);
	    sumb = sosqm(at,N+qd-1);

	    hmax = array_max_abs(h,k);
	    //printf("HMAX %g %g %g \n",hmax,suminc,sumb);

	    if (suminc < sumb) {
	    	if (hmax < e1) {
	    		retval = 1;
	    		break;
	    	} else {
	    		for(i = 0; i < k;++i) {
	    			b[i] = binc[i];
	    		}
	    		lam = lam / vmul;
	    		vmul = 2.0;
				if (M == 1) {
					for(j = 0; j < N;++j) {
						orig[j] = origval[j] - b[0];
					}
				} else {
					for(j = 0; j < N;++j) {
						orig[j] = inp[j];
					}				
				}
			    avaluem(orig,N,p,q,phiinc,thetainc,qd,at);
			    for(i = 0; i < N+qd-1;++i) {
			    	resat[i] = at[i];
			    }
				
				for(i = 0; i < k;++i) {
					stepsize = b[i]*di;//di;//b[i] * di;
					temp = b[i];
					b[i] += stepsize;
					stepsize = b[i] - temp;
		
					if (M == 1) {
						for(j = 0;j < p;++j) {
							phi0[j] = b[1+j];
						}
						for(j = 0; j < q;++j) {
							theta0[j] = b[p+1+j];
						}

						for(j = 0; j < N;++j) {
							orig[j] = origval[j] - b[0];
						}

					} else {
						for(j = 0;j < p;++j) {
							phi0[j] = b[j];
						}
						for(j = 0; j < q;++j) {
							theta0[j] = b[p+j];
						}
						for(j = 0; j < N;++j) {
							orig[j] = inp[j];
						}
					}
				    avaluem(orig,N,p,q,phi0,theta0,qd,at0);	
					t = i * (N + qd - 1);
		
					for(j = 0; j < N + qd - 1;++j) {
						x[t+j] = (at[j] - at0[j]) / stepsize;
					}
		
					b[i] -= stepsize;
				}
				for(i = 0; i < k;++i) {
					t = i * (N + qd - 1);
					for(j = 0; j < k;++j) {
						t2 = j * (N + qd - 1);
						A[i*k+j] = 0.0;
						for(l = 0; l < N + qd -1;++l) {
							A[i*k+j] += x[t+l] * x[t2+l];
						}

					}
				}

				for(i = 0; i < k;++i) {
					t = i * (N + qd - 1);
					g[i] = 0.0;
					for(l = 0; l < N + qd - 1;++l) {
						g[i] += x[t+l] * at[l];
					}
					D[i] = sqrt(A[i*k+i]);
				}
				
				for(i = 0; i < k;++i) {
					for(j = 0; j < k;++j) {
						A[i*k+j] = A[i*k+j] / (D[i] * D[j]);
						if (i == j) {
							A[i*k+j] = 1.0 + lam;
						}
					}
				}

				for(i = 0; i < k; ++i) {
					g[i] = g[i] / D[i];
				}

	    	}
	    } else {
			lam=lam*vmul;
			vmul=2.0*vmul;
	    }

	}
	
	//printf("Iterations %d \n",it);
	
	if (it >= itmax) {
		retval = 4;
	} else if (it < itmax && it > 0) {
		retval = 1;
	}
	// Residual Variance
	*var = sumb * 1.0 /(N - p - q - M);
	//Constant Term
	if (M == 1) {
		temp = 0.0;
		for (i = 0;i < p;++i) {
			temp += phiinc[i];
		}
		*thetac = b[0] * (1.0 - temp);
	} else {
		*thetac = 0.0;
	}
	// Parameters Phi and Theta
	
	for (i = 0; i < p;++i) {
		phi[i] = phiinc[i];
	}
	for (i = 0; i < q;++i) {
		theta[i] = thetainc[i];
	}

	for(i = 0; i < k;++i) {
		t = i * (N + qd - 1);
		for(j = 0; j < k;++j) {
			t2 = j * (N + qd - 1);
			A[i*k+j] = 0.0;
			for(l = 0; l < N + qd -1;++l) {
				A[i*k+j] += x[t+l] * x[t2+l];
			}

		}
	}

	ludecomp(A,k,ipiv);
	minverse(A,k,ipiv,varcovar);

	//mdisplay(varcovar, 1, k*k);

	for(i = 0; i < k*k;++i) {
		varcovar[i] = varcovar[i] * (*var);
	}

	for(i = qd-1; i < N+qd-1;++i) {
		residuals[i-qd+1] = resat[i];
	}


	free(inp);
	free(orig);
	free(phi0);
	free(theta0);
	free(b);
	free(at);
	free(at0);
	free(x);
	free(A);
	free(g);
	free(D);
	free(h);
	free(binc);
	free(ipiv);
	free(phiinc);
	free(thetainc);
	free(atinc);
	free(resat);
	free(origval);
	free(Atemp);
	return retval;
}

int nlalsms(double *vec,int N,int p,int delta, int q, double *phi,double *theta,
		int s,int lps,int deltas,int lqs,double *phis,double *thetas,
		int M,double *thetac,double *var,double eps,double *varcovar,double *residuals) {

	int retval,i,j,l,k,t,t2,it,Norig,qd,itmax;
	double lam,vmul,e1,di,eps2,stepsize,temp,tau,hmax,suminc,sumb,temp2;
	double *inp,*b,*at,*at0,*x,*phi0,*theta0,*orig,*A,*g,*D,*h,*binc,*Atemp;
	double *phiincx,*thetaincx,*atinc,*resat,*origval;
	int *ipiv;
	double *inp2;
	double *phit2,*thetat2,*phist1,*thetast1,*phix,*thetax;
	double *phiinc, *thetainc, *phisinc, *thetasinc;


	eps2 = sqrt(eps);
	lam = 0.01;
	vmul = 2.0;
	e1 = 1.0e-08;
	it = 0;
	Norig = N;
	k = p + q + lps + lqs;
	di = eps2;//0.001;//eps2;//0.001
	qd = q + s*lqs + 20 * (p + s*lps);//15 + deltas*s;
	hmax = 1;
	itmax = 200;
	retval = 0;

	if (qd < (15 + deltas*s)) {
		qd = 15 + deltas*s;
	}

	inp2 = (double*) malloc(sizeof(double) * (N - deltas*s));
	origval = (double*) malloc(sizeof(double) * (N - deltas*s));

	if (deltas > 0) {
		N = diffs(vec,N,deltas,s,inp2);
		for (i = 0; i < N;++i) {
			origval[i] = inp2[i];
		}
	} else {
		for(i = 0; i < N;++i) {
			inp2[i] = vec[i];
			origval[i] = inp2[i];
		}
	}

	inp = (double*) malloc(sizeof(double) * (N - delta));
	orig = (double*)malloc(sizeof(double)* (N - delta));


	if (delta > 0) {
		//free(origval);
		origval = (double*) realloc(origval,sizeof(double) * (N - delta));
		N = diff(inp2,N,delta,inp);
		for (i = 0; i < N;++i) {
			origval[i] = inp[i];
		}
	} else {
		for(i = 0; i < N;++i) {
			inp[i] = inp2[i];
			origval[i] = inp[i];
		}
	}

	USPE(inp,N,p,q,phi,theta,thetac,var); // Initial Parameters
	USPE_seasonal(inp,N,s,lps,lqs,phis,thetas);
	//mdisplay(phi,1,p);
	//mdisplay(theta,1,q);
	//mdisplay(phis,1,lps);
	//mdisplay(thetas,1,lqs);

	// ****Replace this with paramsarima

	// STARTS HERE
/*
	*var = 1.0;

	for(i = 0; i < p;++i) {
		phi[i] = 0.5;
	}

	for(i = 0; i < lps;++i) {
		phis[i] = 0.5;
	}

	for(i = 0; i < q;++i) {
		theta[i] = 0.5;
	}

	for(i = 0; i < lqs;++i) {
		thetas[i] = 0.5;
	}
*/
	//ENDS HERE

	if (M == 1) {
		for(i = 0; i < N;++i) {
			inp[i] -= *thetac;
		}
		M = 1;
		k++;
	} else {
		M = 0;
	}

	b = (double*) malloc(sizeof(double) * k);
	at = (double*) malloc(sizeof(double) * (N+qd-1));
	at0 = (double*) malloc(sizeof(double) * (N+qd-1));
	atinc = (double*) malloc(sizeof(double) * (N+qd-1));
	resat = (double*) malloc(sizeof(double) * (N+qd-1));
	x = (double*) malloc(sizeof(double) * k * (N+qd-1));
	A = (double*) malloc(sizeof(double) * k * k);
	Atemp = (double*) malloc(sizeof(double) * k * k);
	g = (double*) malloc(sizeof(double) * k);
	D = (double*) malloc(sizeof(double) * k);
	h = (double*) malloc(sizeof(double) * k);
	binc = (double*) malloc(sizeof(double) * k);
	ipiv = (int*) malloc(sizeof(int) * k);

	if (M == 1) {
		b[0] = *thetac;
	}

	for(i = 0; i < p;++i) {
		b[M+i] = phi[i];
	}

	for(i = 0; i < q;++i) {
		b[M+p+i] = theta[i];
	}

	for(i = 0; i < lps;++i) {
		b[M+p+q+i] = phis[i];
	}

	for(i = 0; i < lqs;++i) {
		b[M+p+q+lps+i] = thetas[i];
	}

	phi0 = (double*) malloc(sizeof(double) * (p+1));
	theta0 = (double*) malloc(sizeof(double) * (q+1));
	phiinc = (double*)malloc(sizeof(double)* (p + 1));
	thetainc = (double*)malloc(sizeof(double)* (q + 1));
	phiincx = (double*)malloc(sizeof(double)* (s*lps + 1 + p));
	thetaincx = (double*)malloc(sizeof(double)* (s*lqs + 1 + q));
	phist1 = (double*) malloc(sizeof(double) * (lps+1));
	thetast1 = (double*) malloc(sizeof(double) * (lqs+1));
	phit2 = (double*) malloc(sizeof(double) * (s*lps+1));
	thetat2 = (double*) malloc(sizeof(double) * (s*lqs+1));
	phix = (double*) malloc(sizeof(double) * (s*lps+1+p));
	thetax = (double*) malloc(sizeof(double) * (s*lqs+1+q));
	phisinc = (double*)malloc(sizeof(double)* (s*lps + 1));
	thetasinc = (double*)malloc(sizeof(double)* (s*lqs + 1));

	phi0[0] = theta0[0] = phist1[0] = thetast1[0] = 1.0;

	for(i = 0; i < p;++i) {
		phi0[1+i] = - phi[i];
	}

	for(i = 0; i < q;++i) {
		theta0[1+i] = - theta[i];
	}

	for(i = 0; i < lps;++i) {
		phist1[1+i] = - phis[i];
	}

	for(i = 0; i < lqs;++i) {
		thetast1[1+i] = - thetas[i];
	}

	upsample(phist1,lps+1,s,phit2);
	upsample(thetast1,lqs+1,s,thetat2);
	poly(phi0,phit2,phix,1+p,s*lps+1);
	poly(theta0,thetat2,thetax,1+q,s*lqs+1);

	//mdisplay(phit2, 1, s*lps + 1);

	//mdisplay(phit2,1,s*lps+1);
	//mdisplay(thetat2, 1, s*lqs + 1);

	// Edit these Values --

	for (i = 1; i < s*lps + p + 1; ++i) {
		phix[i] = -1.0 * phix[i];
	}

	for (i = 1; i < s*lqs + q + 1; ++i) {
		thetax[i] = -1.0 * thetax[i];
	}

    avaluem(inp,N,s*lps+p,s*lqs+q,phix+1,thetax+1,qd,at);

    //mdisplay(at,1,N+qd-1);

	for (i = 0; i < N + qd - 1; ++i) {
		resat[i] = at[i];
	}
	
	for (i = 0; i < k; ++i) {
		stepsize = b[i] * di;//di;//b[i] * di;
		temp = b[i];
		b[i] += stepsize;
		stepsize = b[i] - temp;
		phi0[0] = theta0[0] = phist1[0] = thetast1[0] = 1.0;
		if (M == 1) {
			for (j = 0; j < p; ++j) {
				phi0[1+j] = - b[1 + j];
			}
			for (j = 0; j < q; ++j) {
				theta0[1+j] = -b[p + 1 + j];
			}

			for (j = 0; j < lps; ++j) {
				phist1[1 + j] = -b[p+q+1+j];
			}

			for (j = 0; j < lqs; ++j) {
				thetast1[1 + j] = -b[p+q+lps+1+j];
			}

			for (j = 0; j < N; ++j) {
				orig[j] = origval[j] - b[0];
			}
		}
		else {
			for (j = 0; j < p; ++j) {
				phi0[1 + j] = -b[j];
			}
			for (j = 0; j < q; ++j) {
				theta0[1 + j] = -b[p + j];
			}

			for (j = 0; j < lps; ++j) {
				phist1[1 + j] = -b[p + q + j];
			}

			for (j = 0; j < lqs; ++j) {
				thetast1[1 + j] = -b[p + q + lps + j];
			}
			for (j = 0; j < N; ++j) {
				orig[j] = inp[j];
			}
		}

		upsample(phist1, lps + 1, s, phit2);
		upsample(thetast1, lqs + 1, s, thetat2);
		poly(phi0, phit2, phix, 1 + p, s*lps + 1);
		poly(theta0, thetat2, thetax, 1 + q, s*lqs + 1);

		//mdisplay(phit2, 1, s*lps + 1);


		for (j = 1; j < s*lps + p + 1; ++j) {
			phix[j] = -1.0 * phix[j];
		}

		for (j = 1; j < s*lqs + q + 1; ++j) {
			thetax[j] = -1.0 * thetax[j];
		}

		avaluem(orig, N, s*lps + p, s*lqs + q, phix + 1, thetax + 1, qd, at0);
		//mdisplay(at0,1,N+qd-1);

		t = i * (N + qd - 1);

		for (j = 0; j < N + qd - 1; ++j) {
			x[t + j] = (at[j] - at0[j]) / stepsize;
		}

		b[i] -= stepsize;
	}

	for (i = 0; i < k; ++i) {
		t = i * (N + qd - 1);
		for (j = 0; j < k; ++j) {
			t2 = j * (N + qd - 1);
			A[i*k + j] = 0.0;
			for (l = 0; l < N + qd - 1; ++l) {
				A[i*k + j] += x[t + l] * x[t2 + l];
			}

		}
	}

	for (i = 0; i < k; ++i) {
		t = i * (N + qd - 1);
		g[i] = 0.0;
		for (l = 0; l < N + qd - 1; ++l) {
			g[i] += x[t + l] * at[l];
		}
		D[i] = sqrt(A[i*k + i]);
	}

	tau = array_max(D, k);

	//printf("TAUMAX %g \n", tau);

	while (hmax >= e1 && it < itmax) {
		it++;
		for (i = 0; i < k; ++i) {
			for (j = 0; j < k; ++j) {
				if (it == 1) {
					A[i*k + j] = A[i*k + j] / (D[i] * D[j]);
				}
				if (i == j) {
					A[i*k + j] = 1.0 + lam;
				}
			}
		}

		if (it == 1) {
			for (i = 0; i < k; ++i) {
				g[i] = g[i] / D[i];
			}
		}

		for(i = 0; i < k*k;++i) {
			Atemp[i] = A[i];
		}
		ludecomp(A,k,ipiv);
		linsolve(A,k,g,ipiv,h);
		for(i = 0; i < k*k;++i) {
			A[i] = Atemp[i];
		}

		for (j = 0; j < k; ++j) {
			h[j] = h[j] / D[j];
			binc[j] = b[j] + h[j];
		}

		phiinc[0] = thetainc[0] = phist1[0] = thetast1[0] = 1.0;

		if (M == 1) {
			for (j = 0; j < p; ++j) {
				phiinc[1 + j] = -binc[1 + j];
			}
			for (j = 0; j < q; ++j) {
				thetainc[1 + j] = -binc[p + 1 + j];
			}

			for (j = 0; j < lps; ++j) {
				phist1[1 + j] = -binc[p + q + 1 + j];
			}

			for (j = 0; j < lqs; ++j) {
				thetast1[1 + j] = -binc[p + q + lps + 1 + j];
			}

			for (j = 0; j < N; ++j) {
				orig[j] = origval[j] - binc[0];
			}
		}
		else {
			for (j = 0; j < p; ++j) {
				phiinc[1 + j] = -binc[j];
			}
			for (j = 0; j < q; ++j) {
				thetainc[1 + j] = -binc[p + j];
			}

			for (j = 0; j < lps; ++j) {
				phist1[1 + j] = -binc[p + q + j];
			}

			for (j = 0; j < lqs; ++j) {
				thetast1[1 + j] = -binc[p + q + lps + j];
			}
			for (j = 0; j < N; ++j) {
				orig[j] = inp[j];
			}
		}

		upsample(phist1, lps + 1, s, phisinc);
		upsample(thetast1, lqs + 1, s, thetasinc);
		poly(phiinc, phisinc, phiincx, 1 + p, s*lps + 1);
		poly(thetainc, thetasinc, thetaincx, 1 + q, s*lqs + 1);

		for (j = 1; j < s*lps + p + 1; ++j) {
			phiincx[j] = -1.0 * phiincx[j];
		}

		for (j = 1; j < s*lqs + q + 1; ++j) {
			thetaincx[j] = -1.0 * thetaincx[j];
		}

		avaluem(orig, N, s*lps + p, s*lqs + q, phiincx + 1, thetaincx + 1, qd, atinc);

		for (i = 0; i < N + qd - 1; ++i) {
			resat[i] = atinc[i];
		}

		suminc = sosqm(atinc, N + qd - 1);
		sumb = sosqm(at, N + qd - 1);

		hmax = array_max_abs(h, k);
		//printf("HMAX %g %g %g \n",hmax,suminc,sumb);

		if (suminc < sumb) {
			if (hmax < e1) {
				retval = 1;
				break;
			}
			else {
				for (i = 0; i < k; ++i) {
					b[i] = binc[i];
				}
				lam = lam / vmul;
				vmul = 2.0;
				if (M == 1) {
					for (j = 0; j < N; ++j) {
						orig[j] = origval[j] - b[0];
					}
				}
				else {
					for (j = 0; j < N; ++j) {
						orig[j] = inp[j];
					}
				}
				avaluem(orig, N, s*lps + p, s*lqs + q, phiincx + 1, thetaincx + 1, qd, at);

				for (i = 0; i < N + qd - 1; ++i) {
					resat[i] = at[i];
				}

				for (i = 0; i < k; ++i) {
					stepsize = b[i] * di;//di;//b[i] * di;
					temp = b[i];
					b[i] += stepsize;
					stepsize = b[i] - temp;

					phi0[0] = theta0[0] = phist1[0] = thetast1[0] = 1.0;
					if (M == 1) {
						for (j = 0; j < p; ++j) {
							phi0[1 + j] = -b[1 + j];
						}
						for (j = 0; j < q; ++j) {
							theta0[1 + j] = -b[p + 1 + j];
						}

						for (j = 0; j < lps; ++j) {
							phist1[1 + j] = -b[p + q + 1 + j];
						}

						for (j = 0; j < lqs; ++j) {
							thetast1[1 + j] = -b[p + q + lps + 1 + j];
						}

						for (j = 0; j < N; ++j) {
							orig[j] = origval[j] - b[0];
						}
					}
					else {
						for (j = 0; j < p; ++j) {
							phi0[1 + j] = -b[j];
						}
						for (j = 0; j < q; ++j) {
							theta0[1 + j] = -b[p + j];
						}

						for (j = 0; j < lps; ++j) {
							phist1[1 + j] = -b[p + q + j];
						}

						for (j = 0; j < lqs; ++j) {
							thetast1[1 + j] = -b[p + q + lps + j];
						}
						for (j = 0; j < N; ++j) {
							orig[j] = inp[j];
						}
					}

					upsample(phist1, lps + 1, s, phit2);
					upsample(thetast1, lqs + 1, s, thetat2);
					poly(phi0, phit2, phix, 1 + p, s*lps + 1);
					poly(theta0, thetat2, thetax, 1 + q, s*lqs + 1);

					//mdisplay(phit2, 1, s*lps + 1);


					for (j = 1; j < s*lps + p + 1; ++j) {
						phix[j] = -1.0 * phix[j];
					}

					for (j = 1; j < s*lqs + q + 1; ++j) {
						thetax[j] = -1.0 * thetax[j];
					}

					avaluem(orig, N, s*lps + p, s*lqs + q, phix + 1, thetax + 1, qd, at0);
					t = i * (N + qd - 1);

					for (j = 0; j < N + qd - 1; ++j) {
						x[t + j] = (at[j] - at0[j]) / stepsize;
					}

					b[i] -= stepsize;
				}
				for (i = 0; i < k; ++i) {
					t = i * (N + qd - 1);
					for (j = 0; j < k; ++j) {
						t2 = j * (N + qd - 1);
						A[i*k + j] = 0.0;
						for (l = 0; l < N + qd - 1; ++l) {
							A[i*k + j] += x[t + l] * x[t2 + l];
						}

					}
				}

				for (i = 0; i < k; ++i) {
					t = i * (N + qd - 1);
					g[i] = 0.0;
					for (l = 0; l < N + qd - 1; ++l) {
						g[i] += x[t + l] * at[l];
					}
					D[i] = sqrt(A[i*k + i]);
				}

				for (i = 0; i < k; ++i) {
					for (j = 0; j < k; ++j) {
						A[i*k + j] = A[i*k + j] / (D[i] * D[j]);
						if (i == j) {
							A[i*k + j] = 1.0 + lam;
						}
					}
				}

				for (i = 0; i < k; ++i) {
					g[i] = g[i] / D[i];
				}

			}
		}
		else {
			lam = lam*vmul;
			vmul = 2.0*vmul;
		}

	}

	if (it >= itmax) {
		retval = 4;
	}
	else if (it < itmax && it > 0) {
		retval = 1;
	}
	// Residual Variance
	*var = sumb * 1.0 / (N - p - q - lps - lqs - M);
	//Constant Term
	if (M == 1) {
		temp = 0.0;
		temp2 = 0.0;
		for (i = 0; i < p; ++i) {
			temp -= phiinc[i+1];
		}
		for (i = 0; i < lps; ++i) {
			temp2 -= phisinc[i + 1];
		}
		*thetac = b[0] * (1.0 - temp)*(1.0 - temp2) ;
	}
	else {
		*thetac = 0.0;
	}
	// Parameters Phi and Theta

	for (i = 0; i < p; ++i) {
		phi[i] = -phiinc[i+1];
	}
	for (i = 0; i < q; ++i) {
		theta[i] = -thetainc[i+1];
	}

	downsample(phisinc, s*lps + 1, s, phist1);
	downsample(thetasinc, s*lqs + 1, s, thetast1);

	//mdisplay(phist1, 1, lps + 1);
	//mdisplay(thetast1, 1, lqs + 1);


	for (i = 0; i < lps; ++i) {
		phis[i] = -phist1[i + 1];
	}

	for (i = 0; i < lqs; ++i) {
		thetas[i] = -thetast1[i + 1];
	}

	for (i = 0; i < k; ++i) {
		t = i * (N + qd - 1);
		for (j = 0; j < k; ++j) {
			t2 = j * (N + qd - 1);
			A[i*k + j] = 0.0;
			for (l = 0; l < N + qd - 1; ++l) {
				A[i*k + j] += x[t + l] * x[t2 + l];
			}

		}
	}

	ludecomp(A, k, ipiv);
	minverse(A, k, ipiv, varcovar);

	for (i = 0; i < k*k; ++i) {
		varcovar[i] = varcovar[i] * (*var);
	}


	for (i = qd - 1; i < N + qd - 1; ++i) {
		residuals[i - qd + 1] = resat[i];
	}
	
	//mdisplay(residuals, 1, N);


	free(inp2);
	free(inp);
	free(orig);
	free(phi0);
	free(theta0);
	free(b);
	free(at);
	free(at0);
	free(x);
	free(A);
	free(g);
	free(D);
	free(h);
	free(binc);
	free(ipiv);
	free(phiincx);
	free(thetaincx);
	free(atinc);
	free(resat);
	free(origval);
	free(phit2);
	free(thetat2);
	free(phist1);
	free(thetast1);
	free(phix);
	free(thetax);
	free(phiinc);
	free(thetainc);
	free(phisinc);
	free(thetasinc);
	free(Atemp);
	return retval;
}
