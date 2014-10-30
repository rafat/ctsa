
#include "emle.h"

alik_css_object alik_css_init(int p, int d, int q, int N) {
	alik_css_object obj = NULL;
	int i, t;
	if (p > 100 || q > 100) {
		printf("\n p and q should be less than 100. \n");
		printf("Use sarima_init if either p > 100 or q > 100 \n");
		exit(1);
	}

	obj = (alik_css_object)malloc(sizeof(struct alik_css_set) + sizeof(double)* 3 * N);

	obj->p = p;
	obj->d = d;
	obj->q = q;
	obj->N = N;
	obj->length = N;
	obj->pq = p + q;
	obj->M = 0;

	if (d == 0) {
		obj->pq = obj->pq + 1;
		obj->M = 1;
	}

	obj->r = p;

	t = 1 + q;
	if (obj->r < t) {
		obj->r = t;
	}
	t = obj->r;
	for (i = 0; i < t; ++i) {
		obj->phi[i] = 0.0;
	}

	for (i = 0; i < t; ++i) {
		obj->theta[i] = 0.0;
	}
	obj->eps = macheps();
	obj->mean = 0.0;

	return obj;
}

alik_object alik_init(int p, int d, int q, int N) {
	alik_object obj = NULL;
	int i,t;
	if (p > 100 || q > 100) {
		printf("\n p and q should be less than 100. \n");
		printf("Use sarima_init if either p > 100 or q > 100 \n");
		exit(1);
	}

	obj = (alik_object)malloc(sizeof(struct alik_set) + sizeof(double) * 3 * N);

	obj->p = p;
	obj->d = d;
	obj->q = q;
	obj->N = N;
	obj->length = N;
	obj->pq = p + q;
	obj->M = 0;

	if (d == 0) {
		obj->pq = obj->pq + 1;
		obj->M = 1;
	}

	obj->r = p;

	t = 1 + q;
	if (obj->r < t) {
		obj->r = t;
	}
	t = obj->r;
	for (i = 0; i < t; ++i) {
		obj->phi[i] = 0.0;
	}

	for (i = 0; i < t; ++i) {
		obj->theta[i] = 0.0;
	}
	obj->eps = macheps();
	obj->mean = 0.0;


	return obj;
}

alik_css_seas_object alik_css_seas_init(int p, int d, int q, int s, int P, int D, int Q, int N) {
	alik_css_seas_object obj = NULL;
	int i, r, t;

	r = p + s*P;

	t = q + s*Q + 1;

	if (r < t) {
		r = t;
	}
	else {
		t = r;
	}

	obj = (alik_css_seas_object)malloc(sizeof(struct alik_css_seas_set) + sizeof(double)* 2 * r + sizeof(double)* 3 * N);

	obj->p = p;
	obj->d = d;
	obj->q = q;
	obj->s = s;
	obj->P = P;
	obj->D = D;
	obj->Q = Q;
	obj->N = N;
	obj->length = N;
	obj->pq = p + q + P + Q;

	obj->M = 0;

	if (d == 0 && D == 0) {
		obj->pq = obj->pq + 1;
		obj->M = 1;
	}

	obj->r = p + s * P;

	t = 1 + q + s*Q;
	if (obj->r < t) {
		obj->r = t;
	}
	t = obj->r;
	obj->offset = 2 * obj->r;

	for (i = 0; i < t; ++i) {
		obj->x[i] = 0.0;// phi
	}

	for (i = 0; i < t; ++i) {
		obj->x[i + r] = 0.0;//theta
	}
	obj->eps = macheps();
	obj->mean = 0.0;

	return obj;
}

alik_seas_object alik_seas_init(int p, int d, int q, int s, int P, int D, int Q, int N) {
	alik_seas_object obj = NULL;
	int i,r, t;

	r = p + s*P;

	t = q + s*Q + 1;

	if (r < t) {
		r = t;
	}
	else {
		t = r;
	}

	obj = (alik_seas_object)malloc(sizeof(struct alik_seas_set) + sizeof(double) * 2 * r + sizeof(double)* 3 * N);

	obj->p = p;
	obj->d = d;
	obj->q = q;
	obj->s = s;
	obj->P = P;
	obj->D = D;
	obj->Q = Q;
	obj->N = N;
	obj->length = N;
	obj->pq = p + q + P + Q;

	obj->M = 0;

	if (d == 0 && D == 0) {
		obj->pq = obj->pq + 1;
		obj->M = 1;
	}

	obj->r = p + s * P;

	t = 1 + q + s*Q;
	if (obj->r < t) {
		obj->r = t;
	}
	t = obj->r;
	obj->offset = 2 * obj->r;

	for (i = 0; i < t; ++i) {
		obj->x[i] = 0.0;// phi
	}

	for (i = 0; i < t; ++i) {
		obj->x[i+r] = 0.0;//theta
	}
	obj->eps = macheps();
	obj->mean = 0.0;

	return obj;
}

static int inclu2(int np, int nrbar, double weight, double *xnext, double *xrow, double ynext,
	double *d, double *rbar, double *thetab,double *ssqerr,double *recres,int *irank) {
	int ifault,i,ithisr,i1,k;
	double y, wt,xi,di,dpi,cbar,sbar,xk,rbthis;

	y = ynext;
	wt = weight;

	for (i = 0; i < np; ++i) {
		xrow[i] = xnext[i];
	}

	*recres = 0.0;
	ifault = 1;
	if (wt <= 0.0) {
		return ifault;
	}
	ifault = 0;
	ithisr = 0;

	for (i = 0; i < np; ++i) {
		if (xrow[i] != 0.0) {
			xi = xrow[i];
			di = d[i];
			dpi = di + wt*xi*xi;
			d[i] = dpi;
			cbar = di / dpi;
			sbar = wt * xi / dpi;
			wt = cbar * wt;
			if (i != (np - 1)) {//40
				i1 = i + 1;
				for (k = i1; k < np; ++k) {
					xk = xrow[k];
					rbthis = rbar[ithisr];
					xrow[k] = xk - xi * rbthis;
					rbar[ithisr] = cbar * rbthis + sbar * xk;
					ithisr++;
				}
			}//40
			xk = y;
			y = xk - xi *thetab[i];
			thetab[i] = cbar * thetab[i] + sbar * xk;
			if (di == 0.0) {
				*irank = *irank + 1;
				return ifault;
			}
		}
		else {
			ithisr += (np - i - 1);
		}
	}
	*ssqerr = *ssqerr + wt*y*y;
	*recres = y * sqrt(wt);

	return ifault;
}

static void regres(int np,int nrbar,double *rbar, double *thetab, double *beta) {
	int ithisr,i,im,i1,jm,j;
	double b1;
	
	ithisr = nrbar-1;
	im = np-1;
	for (i = 0; i < np; ++i) {
		b1 = thetab[im];
		if (im != (np - 1)) {//30
			i1 = i - 1;
			jm = np - 1;
			for (j = 0; j <= i1; ++j) {
				b1 = b1 - rbar[ithisr] * beta[jm];
				ithisr--;
				jm--;
			}
		}//30
		beta[im] = b1;
		im--;
	}
}


int starma(int ip, int iq,double *phi,double *theta,double *A,double *P,double *V) {
	int ifault,ir,np,nrbar,i,ind,j,ir1,irank,ifail;
	int npr, npr1, ind1, ind2, indj,indi,indn;
	double vj,ssqerr,phij,phii,ynext,weight,recres;
	double *thetab, *xnext, *xrow, *rbar;
	ifault = 0;
	/*
	c   algorithm as 154 appl. statist. (1980) vol.29, no.3
	c
	c   invoking this subroutine sets the values of v and phi, and obtains
	c   the initial values of a and p.
	c   this routine is not suitable for use with an ar(1) process.
	c   in this case the following instructions should be used for
	c   initialisation.
	c          v(1)=1.0
	c          a(1)=0.0
	c          p(1)=1.0/(1.0-phi(1)*phi(1))
	c
	*/

	// Check for input errors

	if (ip < 0) {
		ifault = 1;
	}
	if (iq < 0) {
		ifault = 2;
	}

	if (ip*ip + iq*iq == 0) {
		ifault = 4;
	}

	ir = iq + 1;

	if (ir < ip) {
		ir = ip;
	}

	np = (ir * (ir + 1)) / 2;
	nrbar = (np * (np - 1)) / 2;

	if (ir == 1) {
		ifault = 8;
	}
	
	if (ifault != 0) {
		return ifault;
	}

	xnext = (double*)malloc(sizeof(double)* np);
	thetab = (double*)malloc(sizeof(double)* np);
	xrow = (double*)malloc(sizeof(double)* np);
	rbar = (double*)malloc(sizeof(double)* nrbar);

	for (i = 1; i < ir;++i) {
		A[i] = 0.0;
		if (i >= ip) {
			phi[i] = 0.0;
		}
		V[i] = 0.0;
		if (i < iq+1) {
			V[i] = theta[i - 1];
		}
	}

	A[0] = 0.0;
	if (ip == 0) {
		phi[0] = 0.0;
	}
	V[0] = 1.0;
	ind = ir;

	for (j = 1; j < ir; ++j) {
		vj = V[j];
		for (i = j; i < ir; ++i) {
			V[ind] = V[i] * vj;
			ind++;
		}
	}

	if (ip != 0) {//300
		/*
		c   now find p(0).

c   the set of equations s*vec(p(0))=vec(v) is solved for vec(p(0)).
c   s is generated row by row in the array xnext.
c   the order of elements in p is changed, so as to bring more leading
c   zeros into the rows of s, hence achieving a reduction of computing
c   time.
		*/
		ir1 = ir - 1;
		irank = 0;
		ifail = 0;
		ssqerr = 0.0;

		for (i = 0; i < nrbar; ++i) {
			rbar[i] = 0.0;
		}
		for (i = 0; i < np; ++i) {
			P[i] = 0.0;
			thetab[i] = 0.0;
			xnext[i] = 0.0;
		}
		ind = 0;
		ind1 = -1;
		npr = np - ir;
		npr1 = npr + 1;
		indj = npr;
		ind2 = npr-1;

		for (j = 0; j < ir; ++j) {//110
			phij = phi[j];
			xnext[indj] = 0.0;
			indj = indj + 1;
			indi = npr1 + j;

			for (i = j; i < ir; ++i) {//110
				ynext = V[ind];
				ind++;
				phii = phi[i];
				if (j != (ir - 1)) {
					xnext[indj] = -phii;
					if (i != (ir - 1)) {
						xnext[indi] -= phij;
						ind1 ++;
						xnext[ind1] = -1.0;
					}
				}//100
				xnext[npr] = -phii*phij;
				ind2++;
				if (ind2 >= np) {
					ind2 = 0;
				}
				xnext[ind2] += 1.0;
				weight = 1.0;
				ifail = inclu2(np, nrbar, weight, xnext, xrow, ynext, P, rbar, thetab, &ssqerr, &recres, &irank);
				//mdisplay(P, 1, np);
				xnext[ind2] = 0.0;
				if (i != (ir - 1)) {
					xnext[indi] = 0.0;
					indi++;
					xnext[ind1] = 0.0;
				}
			}//110
		}//110
		//mdisplay(P, 1, np);
		regres(np, nrbar, rbar, thetab, P);
		/*
		Now Re-order P
		*/
		ind = npr;
		for (i = 0; i < ir; ++i) {
			ind++;
			xnext[i] = P[ind - 1];
		}
		ind = np;
		ind1 = npr;

		for (i = 0; i < npr; ++i) {
			P[ind - 1] = P[ind1 - 1];
			ind--;
			ind1--;
		}

		for (i = 0; i < ir; ++i) {
			P[i] = xnext[i];
		}
		//return ifault;
	}
	else {

		// P[0] is obtained by back-substitution for a Moving Average process

		indn = np;
		ind = np;
		for (i = 0; i < ir; i++) {
			for (j = 0; j <= i; j++) {
				ind--;
				P[ind] = V[ind];
				if (j != 0)  {
					indn--;
					P[ind] += P[indn];
				}
			}
		}
	}

	free(xnext);
	free(thetab);
	free(xrow);
	free(rbar);
	return ifault;
}


void karma(int ip,int iq,double *phi,double *theta,double *A,double *P,double*V,int N,
	double *W,double *resid,double *sumlog,double *ssq,int iupd,double delta,int *iter,int *nit) {
	int ir,np,i,j,ir1,inde,swtch,ind,indn,l,ii,indw;
	double wnext,dt,A1,ft,ut,g,et;
	double *E;

	ir = iq + 1;
	swtch = 0;
	*iter = 0;

	if (ir < ip) {
		ir = ip;
	}
	E = (double*)malloc(sizeof(double)* ir);
	np = (ir * (ir + 1)) / 2;

	ir1 = ir - 1;
	for (i = 0; i < ir; ++i) {
		E[i] = 0.0;
	}
	inde = 0;

	if (*nit == 0) {
		for (i = 0; i < N; ++i) {//500
			wnext = W[i];
			if (iupd != 1 || i > 0) {//300
				dt = 0.0;
				if (ir != 1) {
					dt = P[ir];
				}
				/*
				if (dt < delta) {
					swtch = 1;
					*nit = i + 1;
					break;//610
				}
				*/
				A1 = A[0];
				if (ir != 1) {//110
					for (j = 0; j < ir1; ++j) {
						A[j] = A[j + 1];
					}
				}//110
				A[ir - 1] = 0.0;
				if (ip != 0) {//200
					for (j = 0; j < ip; ++j) {
						A[j] += phi[j] * A1;
					}
				}//200
				ind = -1;
				indn = ir-1;
				for (l = 0; l < ir; ++l) {
					for (j = l; j < ir; ++j) {
						ind++;
						P[ind] = V[ind];
						if (j != (ir - 1)) {
							indn++;
							P[ind] += P[indn];
						}
					}
				}
			}//300
			ft = fabs(P[0]);//modification
			ut = wnext - A[0];
			if (ir != 1) {//410
				ind = ir;
				for (j = 1; j < ir; ++j) {
					g = P[j] / ft;
					A[j] += g*ut;
					for (l = j; l < ir; ++l) {
						P[ind] -= g*P[l];
						ind++;
					}
				}
			}//410
			A[0] = wnext;
			for (l = 0; l < ir; ++l) {
				P[l] = 0.0;
			}
			resid[i] = ut / sqrt(ft);
			E[inde] = resid[i];
			inde++;
			if (inde >= iq) {
				inde = 0;
			}
			*ssq += (ut*ut) / ft;
			*sumlog += log(ft);// alog
			*iter = *iter + 1;;
			//printf(" sl %g", *sumlog);


		}//500
		if (swtch = 0) {
			*nit = N;
		}
	}
	else {
		i = 0;
		*nit = i;
		for (ii = i; ii < N; ++ii) {//650
			et = W[ii];
			indw = ii;
			if (ip != 0) {//630
				for (j = 0; j < ip; ++j) {
					indw--;
					if (indw >= 0) {
						et -= phi[j] * W[indw];
					}
				}
			}//630
			if (iq != 0) {//645
				for (j = 0; j < iq; ++i) {
					inde--;
					if (inde = -1) {
						inde = iq-1;
					}
					et -= theta[j] * E[inde];
				}
			}//645
			E[inde] = et;
			resid[ii] = et;
			*ssq += et*et;
			*iter = *iter + 1;
			inde++;
			if (inde >= iq) {
				inde = 0;
			}
		}//650

	}
	if (*iter == 0) {
		*iter = 1;
	}
	free(E);
}

int forkal(int ip,int iq,int id,double *phi,double*theta,double *delta,int N,double *W,double *resid,int il,double *Y,double *AMSE) {
	int ifault,ir,np,k,nrbar,ird,irz;
	double *A,*P,*V,*store,*xrow;
	double zero, one, two,AA,del,sumlog,ssq,sigma,A1,dt,phij,phijdt,phii,AMS;
	int i, iq1, j, ll, lli, nt, nj, idk, iid,iupd,nit,iter,ind,irj;
	int ir2, ir1, id2r, id2r1, id1,idd1,idd2,i45,idrr1,iddr,jkl,jkl1,id2r2,ibc,l,iri1,jj;
	int lk, lk1,jklj,iri,kk1,k1,kk,kkk,ind1,ind2,jrj,j1,jrk;
	/*
C     C Translation of Fortran routine
C     ALGORITHM AS 182  APPL. STATIST. (1982) VOL.31, NO.2
C
C     Finite sample prediction from ARIMA processes.
C
C     Auxiliary routines required: KARMA & STARMA from AS 154 and
C     routines called by them: INCLU2 from ASR 17 (a slight variant on
C     AS 75, and REGRES from AS 75.
C*/
	ir = iq + 1;

	if (ir < ip) {
		ir = ip;
	}
	np = (ir * (ir + 1)) / 2;
	nrbar = (np * (np - 1)) / 2;
	ird = ir + id;
	irz = (ird * (ird + 1)) / 2;
	zero = 0.0;
	one = 1.0;
	two = 2.0;

	ifault = 0;

	if (ip < 0) {
		ifault = 1;
	}
	if (iq < 0) {
		ifault = ifault + 2;
	}
	if (ip*ip + iq*iq == 0) {
		ifault = 4;
	}
	if (id < 0) {
		ifault = 8;
	}
	if (il < 1) {
		ifault = 11;
	}
	if (ifault != 0) {
		return ifault;
	}

	A = (double*)malloc(sizeof(double)* ird);
	P = (double*)malloc(sizeof(double)* irz);
	V = (double*)malloc(sizeof(double)* np); 
	store = (double*)malloc(sizeof(double)* ird);
	xrow = (double*)malloc(sizeof(double)* np);

	//Initial Conditions for Kalman Filter

	for (i = 0; i < ird; ++i) {
		A[i] = zero;
		store[i] = zero;
	}
	for (i = 0; i < irz; ++i) {
		P[i] = zero;
	}
	for (i = 0; i < np; ++i) {
		V[i] = zero;
		xrow[i] = zero;
	}
	A[0] = zero;
	V[0] = one;
	/*

	if (np != 1) {
		for (i = 2; i <= np; ++i) {
			V[i-1] = zero;
		}
		if (iq != 0) {
			iq1 = iq + 1;
			for (i = 2; i <= iq1; ++i) {
				V[i-1] = theta[i - 2];
			}
			for (j = 1; j <= iq; ++j) {
				ll = j * (2 * ir + 1 - j) / 2;
				for (i = j; i <= iq; ++i) {
					lli = ll + i;
					V[lli-1] = theta[i-1] * theta[j-1]; // Error
				}
			}
		}
	}//130
*/
	//Find initial likelihood conditions.

	if (ir == 1) {
		*P = 1.0 / (1.0 - phi[0] * phi[0]);
	}
	else {
		starma(ip, iq, phi, theta,A, P, V);
	}

	//Calculate Data Transformations

	nt = N - id;

	if (id != 0) {
		for (j = 1; j <= id; ++j) {
			nj = N - j;
			store[j-1] = W[nj-1];
		}
		for (i = 1; i <= nt; ++i) {
			AA = zero;
			for (k = 1; k <= id; ++k) {
				idk = id + i - k;
				AA -= delta[k - 1] * W[idk - 1];
			}
			iid = i + id;
			W[i - 1] = W[iid - 1] + AA;
		}
	}//170

	//Evaluate likelihood to obtain final KF conditions
	sumlog = ssq = zero;
	del = -1.0;
	iupd = 1;
	nit = 0;
	iter = 0;

	karma(ip, iq, phi, theta, A, P, V, nt, W, resid, &sumlog, &ssq, iupd, del,&iter, &nit);

	//Calculate M.L.E. of sigma squared

	sigma = zero;

	for (j = 0; j < nt; ++j) {
		sigma += resid[j] * resid[j];
	}
	sigma = sigma / nt;
	//mdisplay(resid, 1, N);

	//Reset the initial A and P when differencing occurs

	if (id != 0) {
		for (i = 1; i <= np; ++i) {
			xrow[i-1] = P[i-1];
		}
		for (i = 1; i <= irz; ++i) {
			P[i-1] = zero;
		}
		ind = 0;
		for (j = 1; j <= ir; ++j) {
			k = (j - 1) * (id + ir + 1) - (j - 1) * j / 2;
			for (i = j; i <= ir; ++i) {
				ind++;
				k++;
				P[k - 1] = xrow[ind-1];
			}
		}

		for (j = 1; j <= id; ++j) {
			irj = ir + j;
			A[irj-1] = store[j-1];
		}

	}//250

	//Set up constants
	ir2 = ir + 1;
	ir1 = ir - 1;
	id1 = id - 1;
	id2r = 2 * ird;
	id2r1 = id2r - 1;
	idd1 = 2 * id + 1;
	idd2 = idd1 + 1;
	i45 = id2r + 1;
	idrr1 = ird + 1;
	iddr = 2 * id + ir;
	jkl = ir * (iddr + 1) / 2;
	jkl1 = jkl + 1;
	id2r2 = id2r + 2;
	ibc = ir * (i45 - ir) / 2;

	for (l = 1; l <= il; ++l) {
		//Predict A
		A1 = A[0];
		if (ir != 1) {
			for (i = 1; i <= ir1; ++i) {
				A[i-1] = A[i];
			}
		}//310
		A[ir-1] = zero;
		if (ip != 0) {
			for (j = 1; j <= ip; ++j) {
				A[j-1] += phi[j-1] * A1;
			}
		}//330
		if (id != 0) {
			for (j = 1; j <= id; ++j) {
				irj = ir + j;
				A1 += delta[j-1] * A[irj-1];
			}//340
			if (id >= 2) {
				for (i = 1; i <= id1; ++i) {
					iri1 = ird - i;
					A[iri1] = A[iri1-1];
				}
			}//360
			A[ir2 - 1] = A1;
		}//360

		//A[ir2-1] = A1;

		//Predict P

		if (id != 0) {
			for (i = 1; i <= id; ++i) {
				store[i-1] = zero;
				for (j = 1; j <= id; ++j) {
					ll = imax(i, j);
					k = imin(i, j);
					jj = jkl + (ll - k) + 1 + (k - 1) * (idd2 - k) / 2;
					store[i-1] += delta[j-1] * P[jj-1];
				}
			}//370

			if (id != 1) {
				for (j = 1; j <= id1; ++j) {
					jj = id - j;
					lk = (jj - 1) * (idd2 - jj) / 2 + jkl;
					lk1 = jj * (idd1 - jj) / 2 + jkl;
					for (i = 1; i <= j; ++i) {
						lk = lk + 1;
						lk1 = lk1 + 1;
						P[lk1-1] = P[lk-1];
					}
				}//380
				for (j = 1; j <= id1; ++j) {
					jklj = jkl1 + j;
					irj = ir + j;
					P[jklj - 1] = store[j - 1] + P[irj-1];
				}
			}//400

			P[jkl1 - 1] = P[0];

			for (i = 1; i <= id; ++i) {
				iri = ir + i;
				P[jkl1 - 1] += delta[i - 1] * (store[i - 1] + two * P[iri - 1]);
			}

			for (i = 1; i <= id; ++i) {
				iri = ir + i;
				store[i - 1] = P[iri - 1];
			}

			for (j = 1; j <= ir; ++j) {
				kk1 = j * (id2r1 - j) / 2 + ir;
				k1 = (j - 1) * (id2r - j) / 2 + ir;
				for (i = 1; i <= id; ++i) {
					kk = kk1 + i;
					k = k1 + i;
					P[k - 1] = phi[j - 1] * store[i - 1];
					if (j != ir) {
						P[k - 1] += P[kk - 1];
					}
				}
			}

			for (j = 1; j <= ir; ++j) {
				store[j - 1] = zero;
				kkk = j * (i45 - j) / 2 - id;
				for (i = 1; i <= id; ++i) {
					kkk++;
					store[j - 1] += delta[i - 1] * P[kkk - 1];
				}
			}//440

			if (id != 1) {
				for (j = 1; j <= ir; ++j) {
					k = j * idrr1 - j * (j + 1) / 2 + 1;
					for (i = 1; i <= id1; ++i) {
						k--;
						P[k - 1] = P[k - 2];
					}
				}
			}//460

			for (j = 1; j <= ir; ++j) {
				k = (j - 1) * (id2r - j) / 2 + ir + 1;
				P[k - 1] = store[j - 1] + phi[j - 1] * P[0];
				if (j < ir) {
					P[k - 1] += P[j];
				}
			}

		}//480

		for (i = 0; i < ir; ++i) {
			store[i] = P[i];
		}

		ind = 0;
		dt = P[0];
		for (j = 1; j <= ir; ++j) {
			phij = phi[j - 1];
			phijdt = phij * dt;
			ind2 = (j - 1) * (id2r2 - j) / 2;
			ind1 = j * (i45 - j) / 2;
			for (i = j; i <= ir; ++i) {
				ind++;
				ind2++;
				phii = phi[i - 1];
				P[ind2 - 1] = V[ind - 1] + phii * phijdt;
				if (j < ir) {
					P[ind2 - 1] += store[j] * phii;
				}
				if (i != ir) {
					ind1++;
					P[ind2 - 1] += store[i] * phij + P[ind1 - 1];
				}//500
			}
		}//500

		//Predict Y

		Y[l - 1] = A[0];
		if (id != 0) {
			for (j = 1; j <= id; ++j) {
				irj = ir + j;
				Y[l - 1] += A[irj - 1] * delta[j - 1];
			}
			//Calculate MSE of Y
		}//520
		AMS = P[0];
		if (id != 0) {
			for (j = 1; j <= id; ++j) {
				jrj = ibc + (j - 1) * (idd2 - j) / 2;
				irj = ir + j;
				AMS += (two * delta[j - 1] * P[irj - 1] + P[jrj] * delta[j - 1] * delta[j - 1]);
			}
			if (id != 1) {
				for (j = 1; j <= id1; ++j) {
					j1 = j + 1;
					jrk = ibc + 1 + (j - 1) * (idd2 - j) / 2;
					for (i = j1; i <= id; ++i) {
						jrk++;
						AMS += two * delta[i-1] * delta[j-1] * P[jrk-1];
					}
				}
			}//550
		}//550
		AMSE[l - 1] = AMS * sigma;

	}//560
	/*
	mdisplay(A, 1, ird);
	mdisplay(P, 1, irz);
	mdisplay(V, 1, np);
	mdisplay(store, 1, ird);
	mdisplay(xrow, 1, np);
	*/
	free(A);
	free(P);
	free(V);
	free(store);
	free(xrow);
	return ifault;
}

double fcss(double *b, int pq, void *params) {
	double value,ssq,temp;
	int ip, iq, N,i,ncond,j,jm,iter;
	double *phi, *theta;
	alik_css_object obj = (alik_css_object)params;

	ip = obj->p;
	iq = obj->q;

	N = obj->N;

	value = ssq = 0.0;
	ncond = ip;
	iter = 0;
	phi = (double*)malloc(sizeof(double)* ip);
	theta = (double*)malloc(sizeof(double)* iq);

	for (i = 0; i < ip; ++i) {
		phi[i] = b[i];
	}

	for (i = 0; i < iq; ++i) {
		theta[i] = b[i + ip];
	}

	for (i = 0; i < ncond; ++i) {
		obj->x[N + i] = 0.0;
	}

	if (obj->M == 1) {
		for (i = 0; i < N; ++i) {
			obj->x[i] = obj->x[2 * N + i] - b[ip+iq];
		}
	}

	for (i = ncond; i < N; ++i) {
		iter++;
		temp = obj->x[i];

		for (j = 0; j < ip; ++j) {
			temp = temp -  phi[j] * obj->x[i - j - 1];
		}

		if (i - ncond < iq) {
			jm = i - ncond;
		}
		else {
			jm = iq;
		}

		for (j = 0; j < jm; ++j) {
			temp = temp - theta[j] * obj->x[N + i - j - 1];
		}

		obj->x[N + i] = temp;
		ssq += temp * temp;
	}
	obj->ssq = ssq;
	value = 0.5 * log(ssq / (double)iter);
	obj->loglik = value;

	free(phi);
	free(theta);
	return value;
}

int css(double *inp, int N, int optmethod, int p, int d, int q, double *phi, double *theta, double *wmean,double *var,double *resid,double *loglik,double *hess) {
	int i, pq, retval, length,M,ret;
	double *b, *tf, *x,*dx,*thess;
	int *ipiv;
	double maxstep;
	alik_css_object obj;

	x = (double*)malloc(sizeof(double)* (N - d));


	length = N;

	maxstep = 1.0;

	obj = alik_css_init(p, d, q, N);
	pq = obj->pq;
	b = (double*)malloc(sizeof(double)* pq);
	tf = (double*)malloc(sizeof(double)* pq);
	thess = (double*)malloc(sizeof(double)* pq*pq);
	dx = (double*)malloc(sizeof(double)* pq);
	ipiv = (int*)malloc(sizeof(int)* pq);
	/*

	*/

	if (d > 0) {
		N = diff(inp, N, d, x); // No need to demean x
		M = 0;
		*wmean = 0.0;
	}
	else {
		*wmean = mean(inp, N);
		for (i = 0; i < N; ++i) {
			x[i] = inp[i];
		}
		M = 1;
	}

	obj->N = N;
	obj->mean = *wmean;


	for (i = 0; i < p; ++i) {
		b[i] = 0.0;
	}
	for (i = 0; i < q; ++i) {
		b[p + i] = 0.0;
	}

	if (obj->M == 1) {
		b[p + q] = obj->mean;
	}

	for (i = 0; i < N; ++i) {
		obj->x[i] = obj->x[2*N+i] = x[i];
	}
	for (i = N; i < 2 * N; ++i) {
		obj->x[i] = 0.0;
	}

	custom_function fcss_min = { fcss, obj };
	retval = fminunc(&fcss_min, NULL, pq, b, maxstep, optmethod, tf);

	if (retval == 0) {
		ret = 0;
	} else if (retval == 15) {
		ret = 15;
	} else if (retval == 4) {
		ret = 4;
	}
	else {
		ret = 1;
	}

	for (i = 0; i < pq; ++i) {
		dx[i] = 1.0;
	}

	hessian_fd(&fcss_min, tf, pq, dx, obj->eps, hess);
	//mdisplay(hess, pq, pq);
	mtranspose(hess, pq, pq, thess);

	for (i = 0; i < pq*pq; ++i) {
		thess[i] = (N - d) * 0.5 * (hess[i] + thess[i]);
	}


	ludecomp(thess, pq, ipiv);
	minverse(thess, pq, ipiv, hess);

	for (i = 0; i < p; ++i) {
		phi[i] = tf[i];
	}
	for (i = 0; i < q; ++i) {
		theta[i] = -tf[p + i];
	}
	if (obj->M == 1) {
		*wmean = tf[p + q];
	}
	else {
		*wmean = 0.0;
	}
	/*
	if (M == 1) {
		temp = 0.0;
		for (i = 0; i < p; ++i) {
			temp += phi[i];
		}
		*wmean = *wmean * (1.0 - temp);
	}
	else {
		*wmean = 0.0;
	}
	*/
	*var = (obj->ssq) / (double)N;
	for (i = 0; i < N - d; ++i) {
		resid[i] = obj->x[N+i];
	}
	*loglik = obj->loglik;
	//printf("MEAN %g \n", mean(obj->x+N,N));
	//mdisplay(obj->x + N, 1, N);

	free(b);
	free(tf);
	free(x);
	free(thess);
	free(dx);
	free(ipiv);
	free_alik_css(obj);
	return ret;
}


double fas154(double *b,int pq,void *params) {
	double value,ssq,sumlog,delta;
	int ip, iq, ir,i,np,ifault,N,iupd,nit,iter;
	double *phi, *theta, *A,*P,*V;

	value = ssq = sumlog = 0.0;
	alik_object obj = (alik_object)params;

	ip = obj->p;
	iq = obj->q;
	ir = obj->r;
	N = obj->N;
	np = (ir * (ir + 1)) / 2;

	phi = (double*)malloc(sizeof(double)* ir);
	theta = (double*)malloc(sizeof(double)* ir);
	A = (double*)malloc(sizeof(double)* ir);
	P = (double*)malloc(sizeof(double)* np);
	V = (double*)malloc(sizeof(double)* np);
	//bt = (double*)malloc(sizeof(double)* (ip+iq));
	/*
	for (i = 0; i < ip + iq; ++i) {
		bt[i] = b[i];
	}

	pdlreg(ip, bt, b);
	pdlreg(iq, bt + ip, b + ip);
	*/
	for (i = 0; i < ip; ++i) {
		phi[i] = b[i];
	}
	for (i = ip; i < ir; ++i) {
		phi[i] = 0.0;
	}

	for (i = 0; i < iq; ++i) {
		theta[i] = b[i+ip];
	}
	for (i = iq; i < ir; ++i) {
		theta[i] = 0.0;
	}
	if (obj->M == 1) {
		for (i = 0; i < N; ++i) {
			obj->x[i] = obj->x[2 * N + i] - b[ip + iq];
		}
	}

	iupd = 0;

	if (ip == 1 && iq == 0) {
		*V = 1.0;
		*A = 0.0;
		*P = 1.0 / (1.0 - phi[0] * phi[0]);
	}
	else {
		iupd = 1;
		ifault = starma(ip, iq, phi, theta, A, P, V);
	}
	
	nit = 0;
	delta = 0.001;
	//mdisplay(P, 1, np);
	karma(ip, iq, phi, theta, A, P, V, N, obj->x, obj->x + N, &sumlog, &ssq, iupd, delta,&iter, &nit);
	obj->ssq = ssq;
	//mdisplay(phi, 1, ip);
	value = 0.5 * (sumlog/(double)iter + log(ssq/(double) iter));
	obj->loglik = value;

	//printf("sumlog ssq %g %g %d \n", sumlog,value,iter);

	free(phi);
	free(theta);
	free(A);
	free(P);
	free(V);
	//free(bt);
	return value;
}

int as154(double *inp, int N, int optmethod, int p, int d, int q, double *phi, double *theta, double *wmean, double *var,double *resid,double *loglik,double *hess) {
	int i,pq,retval,length,ret;
	double *b,*tf,*x,*dx,*thess;
	int *ipiv;
	double maxstep;
	alik_object obj;
	//custom_function as154_min;

	x = (double*)malloc(sizeof(double)* (N - d));

	length = N;
	
	maxstep = 1.0;

	obj = alik_init(p, d, q, N);


	css(inp, N, optmethod, p, d, q, phi, theta, wmean, var,resid,loglik,hess);

	if (d > 0) {
		N = diff(inp, N, d, x); // No need to demean x
	}
	else {
		for (i = 0; i < N; ++i) {
			x[i] = inp[i];
		}
	}

	obj->N = N;
	obj->mean = *wmean;
	pq = obj->pq;
	b = (double*)malloc(sizeof(double)* pq);
	tf = (double*)malloc(sizeof(double)* pq);
	thess = (double*)malloc(sizeof(double)* pq*pq);
	dx = (double*)malloc(sizeof(double)* pq);
	ipiv = (int*)malloc(sizeof(int)* pq);

	for (i = 0; i < p; ++i) {
		b[i] = phi[i];
	}
	for (i = 0; i < q; ++i) {
		b[p + i] = -theta[i];
	}

	if (obj->M == 1) {
		b[p + q] = obj->mean;
	}


	for (i = 0; i < N; ++i) {
		obj->x[i] = obj->x[2 * N + i] = x[i];
	}
	for (i = N; i < 2 * N; ++i) {
		obj->x[i] = 0.0;
	}

	custom_function as154_min = { fas154, obj };
	retval = fminunc(&as154_min, NULL, pq, b,maxstep, optmethod, tf);

	if (retval == 0) {
		ret = 0;
	}
	else if (retval == 15) {
		ret = 15;
	}
	else if (retval == 4) {
		ret = 4;
	}
	else {
		ret = 1;
	}

	for (i = 0; i < pq; ++i) {
		dx[i] = 1.0;
	}

	hessian_fd(&as154_min, tf, pq, dx, obj->eps, hess);
	
	mtranspose(hess, pq, pq, thess);

	for (i = 0; i < pq*pq; ++i) {
		thess[i] = (N-d) * 0.5 * (hess[i] + thess[i]);
	}
	

	ludecomp(thess, pq, ipiv);
	minverse(thess, pq, ipiv, hess);


	for (i = 0; i < p; ++i) {
		phi[i] = tf[i];
	}
	for (i = 0; i < q; ++i) {
		theta[i] = -tf[p + i];
	}
	if (obj->M == 1) {
		*wmean = tf[p + q];
	}
	else {
		*wmean = 0.0;
	}
	/*
	wmean = 0.0;
	for (i = 0; i < N; ++i) {
		wmean += (obj->x[N + i] * obj->x[N + i]);
	}*/

	*var = (obj->ssq) / (double) N;
	for (i = 0; i < N - d; ++i) {
		resid[i] = obj->x[N + i];
	}
	*loglik = obj->loglik;
	//printf("MEAN %g \n", mean(obj->x+N,N));
	//mdisplay(obj->x + N, 1, N);

	free(b);
	free(tf);
	free(x);
	free(thess);
	free(ipiv);
	free(dx);
	free_alik(obj);
	return ret;
}

double fcss_seas(double *b, int pq, void *params) {
	double value, ssq,temp;
	int p, q, ps, qs, s, offset,jm;
	int ip, iq, i, j, N,iter,ncond;
	double *phi, *theta;

	value = ssq =  0.0;
	
	alik_css_seas_object obj = (alik_css_seas_object)params;

	ip = obj->p + (obj->s * obj->P);
	iq = obj->q + (obj->s * obj->Q);

	N = obj->N;
	p = obj->p;
	ps = obj->P;
	q = obj->q;
	qs = obj->Q;
	s = obj->s;
	offset = obj->offset;
	ncond = p + s* ps;

	//printf("\n %d %d %d \n", s,ps, qs);

	phi = (double*)malloc(sizeof(double)* ip);
	theta = (double*)malloc(sizeof(double)* iq);



	for (i = 0; i < p; ++i) {
		phi[i] = b[i];
	}

	for (i = p; i < ip; ++i) {
		phi[i] = 0.0;
	}

	for (i = 0; i < q; ++i) {
		theta[i] = b[i + p];
	}

	for (i = q; i < iq; ++i) {
		theta[i] = 0.0;
	}

	for (j = 0; j < ps; ++j) {
		phi[(j + 1)*s - 1] += b[p + q + j];
		for (i = 0; i < p; ++i) {
			phi[(j + 1)*s + i] -= b[i] * b[p + q + j];
		}
	}

	for (j = 0; j < qs; ++j) {
		theta[(j + 1)*s - 1] += b[p + q + ps + j];
		for (i = 0; i < q; ++i) {
			theta[(j + 1)*s + i] += b[i + p] * b[p + q + ps + j];
		}
	}
	//mdisplay(theta, 1, q + s*qs);
	//mdisplay(b, 1, p + q + ps + qs);
	for (i = 0; i < ncond; ++i) {
		obj->x[offset+N + i] = 0.0;
	}

	if (obj->M == 1) {
		for (i = 0; i < N; ++i) {
			obj->x[offset+i] = obj->x[offset+ 2 * N + i] - b[p+q+ps+qs];
		}
	}
	iter = 0;
	for (i = ncond; i < N; ++i) {
		iter++;
		temp = obj->x[offset+i];

		for (j = 0; j < ip; ++j) {
			temp = temp - phi[j] * obj->x[offset + i - j - 1];
		}

		if (i - ncond < iq) {
			jm = i - ncond;
		}
		else {
			jm = iq;
		}

		for (j = 0; j < jm; ++j) {
			temp = temp - theta[j] * obj->x[offset+N + i - j - 1];
		}

		obj->x[offset+N + i] = temp;
		ssq += temp * temp;
	}
	obj->ssq = ssq;
	value = 0.5 * log(ssq / (double)iter);
	obj->loglik = value;

	//printf("ssq %g %g \n", value,ssq);

	free(phi);
	free(theta);
	return value;
}

int css_seas(double *inp, int N, int optmethod, int p, int d, int q, int s, int P, int D, int Q,
	double *phi, double *theta, double *PHI, double *THETA, double *wmean, double *var,double *loglik,double *hess) {
	int i, pq, retval, length, offset,ret;
	double *b, *tf, *x, *inp2,*dx,*thess;
	int *ipiv;
	double maxstep;
	alik_css_seas_object obj;
	//custom_function as154_min;

	inp2 = (double*)malloc(sizeof(double)* (N - s*D));

	length = N;

	maxstep = 1.0;

	obj = alik_css_seas_init(p, d, q, s, P, D, Q, N);

	pq = obj->pq;
	b = (double*)malloc(sizeof(double)* pq);
	tf = (double*)malloc(sizeof(double)* pq);
	thess = (double*)malloc(sizeof(double)* pq*pq);
	dx = (double*)malloc(sizeof(double)* pq);
	ipiv = (int*)malloc(sizeof(int)* pq);
	/*

	*/

	if (D > 0) {
		N = diffs(inp, N, D, s, inp2);
	}
	else {
		for (i = 0; i < N; ++i) {
			inp2[i] = inp[i];
		}
	}

	x = (double*)malloc(sizeof(double)* (N - d));

	if (d > 0) {
		N = diff(inp2, N, d, x); // No need to demean x
	}
	else {
		for (i = 0; i < N; ++i) {
			x[i] = inp2[i];
		}
	}

	obj->N = N;

	offset = obj->offset;
	// Test

	for (i = 0; i < p; ++i) {
		b[i] = 0.0;
	}
	for (i = 0; i < q; ++i) {
		b[p + i] = 0.0;
	}
	for (i = 0; i < P; ++i) {
		b[p + q + i] = 0.0;
	}
	for (i = 0; i < Q; ++i) {
		b[p + q + P + i] = 0.0;
	}

	if (obj->M == 1) {
		b[p + q + P + Q] = mean(x, N);
		*wmean = b[p + q + P + Q];
	}
	else {
		*wmean = 0.0;
	}

	obj->mean = *wmean;

	//mdisplay(b, 1, p + q + P + Q);

	for (i = 0; i < N; ++i) {
		obj->x[offset + i] = obj->x[offset + 2 * N + i] = x[i];
	}
	for (i = N; i < 2 * N; ++i) {
		obj->x[offset + i] = 0.0;
	}

	custom_function css_min = { fcss_seas, obj };
	retval = fminunc(&css_min, NULL, pq, b, maxstep, optmethod, tf);
	if (retval == 0) {
		ret = 0;
	}
	else if (retval == 15) {
		ret = 15;
	}
	else if (retval == 4) {
		ret = 4;
	}
	else {
		ret = 1;
	}

	for (i = 0; i < pq; ++i) {
		dx[i] = 1.0;
	}

	hessian_fd(&css_min, tf, pq, dx, obj->eps, hess);

	mtranspose(hess, pq, pq, thess);

	for (i = 0; i < pq*pq; ++i) {
		thess[i] = (N - d - s*D) * 0.5 * (hess[i] + thess[i]);
	}


	ludecomp(thess, pq, ipiv);
	minverse(thess, pq, ipiv, hess);

	for (i = 0; i < p; ++i) {
		phi[i] = tf[i];
	}
	for (i = 0; i < q; ++i) {
		theta[i] = -tf[p + i];
	}
	for (i = 0; i < P; ++i) {
		PHI[i] = tf[p + q + i];
	}
	for (i = 0; i < Q; ++i) {
		THETA[i] = -tf[p + q + P + i];
	}

	if (obj->M == 1) {
		*wmean = tf[p + q + Q + P];
	}
	else {
		*wmean = 0.0;
	}
	
	/*
	wmean = 0.0;
	for (i = 0; i < N; ++i) {
	wmean += (obj->x[N + i] * obj->x[N + i]);
	}*/

	*var = (obj->ssq) / (double)N;
	*loglik = obj->loglik;
	//printf("MEAN %g \n", mean(obj->x+N,N));
	//mdisplay(obj->x + N, 1, N);

	free(b);
	free(tf);
	free(inp2);
	free(x);
	free(thess);
	free(dx);
	free(ipiv);
	free_alik_css_seas(obj);
	return ret;
}

double fas154_seas(double *b, int pq, void *params) {
	double value, ssq, sumlog, delta;
	int p, q, ps, qs,s,offset;
	int ip, iq, ir, i,j, np, ifault, N, iupd, nit, iter;
	double *phi, *theta, *A, *P, *V;

	value = ssq = sumlog = 0.0;
	alik_seas_object obj = (alik_seas_object)params;

	ip = obj->p + (obj->s * obj->P);
	iq = obj->q + (obj->s * obj->Q);
	ir = obj->r;
	N = obj->N;
	p = obj->p;
	ps = obj->P;
	q = obj->q;
	qs = obj->Q;
	s = obj->s;
	offset = obj->offset;
	np = (ir * (ir + 1)) / 2;

	phi = (double*)malloc(sizeof(double)* ir);
	theta = (double*)malloc(sizeof(double)* ir);
	A = (double*)malloc(sizeof(double)* ir);
	P = (double*)malloc(sizeof(double)* np);
	V = (double*)malloc(sizeof(double)* np);

	for (i = 0; i < ir; ++i) {
		phi[i] = 0.0;
		theta[i] = 0.0;
	}

	for (i = 0; i < p; ++i) {
		phi[i] = b[i];
	}

	for (i = 0; i < q; ++i) {
		theta[i] = b[i + p];
	}
	/*
	for (i = p; i < ip; ++i) {
		phi[i] = 0.0;
	}

	for (i = q; i < iq; ++i) {
		theta[i] = 0.0;
	}
	*/

	for (j = 0; j < ps; ++j) {
		phi[(j + 1)*s - 1] += b[p + q + j];
		for (i = 0; i < p; ++i) {
			phi[(j + 1)*s + i] -= b[i] * b[p + q + j];
		}
	}

	for (j = 0; j < qs; ++j) {
		theta[(j + 1)*s - 1] += b[p + q + ps + j];
		for (i = 0; i < q; ++i) {
			theta[(j + 1)*s + i] += b[i + p] * b[p + q + ps + j];
		}
	}
	if (obj->M == 1) {
		for (i = 0; i < N; ++i) {
			obj->x[offset + i] = obj->x[offset + 2 * N + i] - b[p + q + ps + qs];
		}
	}
	//mdisplay(phi, 1, ir);
	//mdisplay(theta, 1, ir);
	iupd = 0;
	//mdisplay(b, 1, pq);

	if (ip == 1 && iq == 0) {
		*V = 1.0;
		*A = 0.0;
		*P = 1.0 / (1.0 - phi[0] * phi[0]);
	}
	else {
		iupd = 1;
		ifault = starma(ip, iq, phi, theta, A, P, V);
	}

	nit = 0;
	delta = 0.001;
	//mdisplay(P, 1, np);
	karma(ip, iq, phi, theta, A, P, V, N, obj->x + offset, obj->x + offset + N, &sumlog, &ssq, iupd, delta, &iter, &nit);
	obj->ssq = ssq;
	//mdisplay(V, 1, np);

	value = 0.5 * (sumlog / (double)iter + log(ssq / (double)iter));
	obj->loglik = value;
	//printf("sumlog ssq %g %g %d \n", sumlog,value,iter);

	free(phi);
	free(theta);
	free(A);
	free(P);
	free(V);
	return value;
}

int as154_seas(double *inp, int N, int optmethod, int p, int d, int q, int s, int P, int D, int Q,
	double *phi, double *theta, double *PHI, double *THETA, double *wmean,double *var,double *loglik,double *hess) {
	int i, pq, retval, length, offset,ret;
	double *b, *tf, *x,*inp2,*dx,*thess;
	int *ipiv;
	double maxstep;
	alik_seas_object obj;
	//custom_function as154_min;
	obj = alik_seas_init(p, d, q, s, P, D, Q, N);
	inp2 = (double*)malloc(sizeof(double)* (N - s*D));
	pq = obj->pq;
	b = (double*)malloc(sizeof(double)* pq);
	tf = (double*)malloc(sizeof(double)* pq);
	thess = (double*)malloc(sizeof(double)* pq*pq);
	dx = (double*)malloc(sizeof(double)* pq);
	ipiv = (int*)malloc(sizeof(int)* pq);

	length = N;

	maxstep = 1.0;


	css_seas(inp, N, optmethod, p, d, q, s, P, D, Q, phi, theta, PHI, THETA, wmean, var,loglik,hess);

	/*

	*/

	if (D > 0) {
		N = diffs(inp, N, D, s, inp2);
	}
	else {
		for (i = 0; i < N; ++i) {
			inp2[i] = inp[i];
		}
	}

	x = (double*)malloc(sizeof(double)* (N - d));

	if (d > 0) {
		N = diff(inp2, N, d, x); // No need to demean x
	}
	else {
		for (i = 0; i < N; ++i) {
			x[i] = inp2[i];
		}
	}

	obj->N = N;

	offset = obj->offset;
	for (i = 0; i < p; ++i) {
		b[i] = phi[i];
	}
	for (i = 0; i < q; ++i) {
		b[p + i] = -theta[i];
	}
	for (i = 0; i < P; ++i) {
		b[p + q + i] = PHI[i];
	}
	for (i = 0; i < Q; ++i) {
		b[p + q + P + i] = -THETA[i];
	}

	if (obj->M == 1) {
		b[p + q + P + Q] = *wmean;
	}

	obj->mean = *wmean;

	//mdisplay(b, 1, p + q + P + Q);

	for (i = 0; i < N; ++i) {
		obj->x[offset + i] = obj->x[offset + 2 * N + i] = x[i];
	}
	for (i = N; i < 2 * N; ++i) {
		obj->x[offset + i] = 0.0;
	}
	//printf("\n %d %g ", pq,maxstep);

	custom_function as154_min = { fas154_seas, obj };
	retval = fminunc(&as154_min, NULL, pq, b, maxstep, optmethod, tf);
	if (retval == 0) {
		ret = 0;
	}
	else if (retval == 15) {
		ret = 15;
	}
	else if (retval == 4) {
		ret = 4;
	}
	else {
		ret = 1;
	}

	for (i = 0; i < pq; ++i) {
		dx[i] = 1.0;
	}

	hessian_fd(&as154_min, tf, pq, dx, obj->eps, hess);

	mtranspose(hess, pq, pq, thess);

	for (i = 0; i < pq*pq; ++i) {
		thess[i] = (N - d - s*D) * 0.5 * (hess[i] + thess[i]);
	}


	ludecomp(thess, pq, ipiv);
	minverse(thess, pq, ipiv, hess);

	for (i = 0; i < p; ++i) {
		phi[i] = tf[i];
	}
	for (i = 0; i < q; ++i) {
		theta[i] = -tf[p + i];
	}
	for (i = 0; i < P; ++i) {
		PHI[i] = tf[p + q + i];
	}
	for (i = 0; i < Q; ++i) {
		THETA[i] = -tf[p + q + P + i];
	}

	if (obj->M == 1) {
		*wmean = tf[p + q + Q + P];
	}
	else {
		*wmean = 0.0;
	}

	*var = (obj->ssq) / (double) N;
	*loglik = obj->loglik;
	//printf("MEAN %g \n", mean(obj->x+N,N));
	//mdisplay(obj->x + N, 1, N);

	free(b);
	free(tf);
	free(inp2);
	free(x);
	free(dx);
	free(thess);
	free(ipiv);
	free_alik_seas(obj);
	return ret;
}

int flikam(double *P, int MP, double *Q, int MQ,double *W, double *E, int N,double *ssq,double *fact,double *VW,double *VL,int MRP1,double *VK,int MR,double TOLER) {
	int ifault,MXPQ,MXPQP1,MQP1,MPP1,FLAG;
	int k,j,jp2mk,jp1mk,LAST,LOOP,JFROM,i,NEXTI,imj;
	double zero, p0625, one, two, four, sixteen,epsil1;
	double A, ALF, AOR, DETCAR, DETMAN, FLJ, R, VL1, VW1;

	zero = 0.0; p0625 = 0.0625; one = 1.0; two = 2.0; four = 4.0; sixteen = 16.0;
	epsil1 = 1.0e-10;

	*ssq = *fact = zero;
	DETMAN = one;
	DETCAR = zero;
	MXPQ = imax(MP, MQ);// Check
	MXPQP1 = MXPQ + 1;
	MQP1 = MQ + 1;
	MPP1 = MP + 1;
	FLAG = 0;

	ifault = twacf(P, MP,Q, MQ,VW,MXPQP1,VL, MXPQP1,VK,MXPQ);
	if (MR != imax(MP, MQP1)) {
		ifault = 6;
	}
	if (MRP1 != MR + 1) {
		ifault = 7;
	}

	if (ifault > 0) {
		return ifault;
	}

	VK[0] = VW[0];

	if (MR != 1) {
		for (k = 2; k <= MR; ++k) {
			VK[k - 1] = zero;
			if (k <= MP) {
				for (j = k; j <= MP; ++j) {
					jp2mk = j + 2 - k;
					VK[k - 1] += P[j - 1] * VW[jp2mk - 1];
				}
			}//120
			if (k <= MQP1) {
				for (j = k; j <= MQP1; ++j) {
					jp1mk = j + 1 - k;
					VK[k - 1] -= Q[j - 2] * VL[jp1mk - 1];
				}
			}
		}
	}//150

	R = VK[0];
	VL[MR - 1] = zero;
	for (j = 0; j < MR; ++j) {
		VW[j] = zero;
		if (j != MR - 1) {
			VL[j] = VK[j + 1];
		}
		if (j <= MP - 1) {
			VL[j] += P[j] * R;
		}
		VK[j] = VL[j];
	}

	LAST = MPP1 - MQ;
	LOOP = MP;
	JFROM = MPP1;
	VW[MPP1-1] = zero;
	VL[MXPQP1-1] = zero;

	if (N <= 0) {
		return 9;
	}

	for (i = 1; i <= N; ++i) {
		if (i == LAST) {
			LOOP = imin(MP, MQ);
			JFROM = LOOP + 1;
			if (MQ <= 0) {
				FLAG = 1;
				break;
			}
		}
		if (R <= epsil1) {
			return 8;
		}
		if (fabs(R - one) < TOLER && i > MXPQ) {
			FLAG = 1;
			break;
		}
		DETMAN *= R;
		while (fabs(DETMAN) >= one) {
			DETMAN *= p0625;
			DETCAR += four;
		}
		while (fabs(DETMAN) < p0625) {
			DETMAN *= sixteen;
			DETCAR -= four;
		}
		VW1 = VW[0];
		A = W[i - 1] - VW1;
		E[i - 1] = A / sqrt(R);
		AOR = A / R;
		*ssq += A*AOR;
		VL1 = VL[0];
		ALF = VL1 / R;
		R -= ALF * VL1;
		if (LOOP != 0) {
			for (j = 0; j < LOOP; ++j) {
				FLJ = VL[j + 1] + P[j] * VL1;
				VW[j] = VW[j + 1] + P[j] * VW1 + AOR * VK[j];
				VL[j] = FLJ - ALF * VK[j];
				VK[j] -= ALF * FLJ;
			}
		}

		if (JFROM <= MQ) {
			for (j = JFROM; j <= MQ; ++j) {
				VW[j - 1] = VW[j] + AOR * VK[j - 1];
				VL[j - 1] = VL[j] - ALF * VK[j - 1];
				VK[j - 1] -= ALF*VL[j];
			}
		}
		if (JFROM <= MP) {
			for (j = JFROM; j <= MP; ++j) {
				VW[j - 1] = VW[j] + P[j-1] * W[i-1];
			}
		}
	}

	if (FLAG == 1) {
		NEXTI = i;
		ifault = -NEXTI;
		for (i = NEXTI; i <= N; ++i) {
			E[i - 1] = W[i - 1];
		}
		if (MP != 0) {
			for (i = NEXTI; i <= N; ++i) {
				for (j = 1; j <= MP; ++j) {
					imj = i - j;
					E[i - 1] -= P[j - 1] * W[imj - 1];
				}
			}
		}//340
		if (MQ != 0) {
			for (i = NEXTI; i <= N; ++i) {
				for (j = 1; j <= MQ; ++j) {
					imj = i - j;
					E[i - 1] += Q[j - 1] * E[imj - 1];
				}
			}
		}
		for (i = NEXTI; i <= N; ++i) {
			*ssq += E[i - 1] * E[i - 1];
		}
	}
	/*
	mdisplay(VW, 1, MRP1);
	mdisplay(VL, 1, MRP1);
	mdisplay(VK, 1, MR);
	*/
	*fact = pow(DETMAN, one / (double)N) * pow(two, DETCAR / (double)N);

	return ifault;
}

double fas197(double *b, int pq, void *params) {
	double value, ssq, fact, delta;
	int ip, iq, i, ifault, N;
	double *phi, *theta,*VW,*VK,*VL;
	int MR, MRP1;

	value = ssq = fact = 0.0;
	alik_object obj = (alik_object)params;

	ip = obj->p;
	iq = obj->q;
	N = obj->N;
	MR = imax(ip, iq + 1);
	MRP1 = MR + 1;
	phi = (double*)malloc(sizeof(double)* ip);
	theta = (double*)malloc(sizeof(double)* iq);
	VW = (double*)malloc(sizeof(double*)*MRP1);
	VL = (double*)malloc(sizeof(double*)*MRP1);
	VK = (double*)malloc(sizeof(double*)*MR);

	for (i = 0; i < ip; ++i) {
		phi[i] = b[i];
	}

	for (i = 0; i < iq; ++i) {
		theta[i] = b[i + ip];
	}

	if (obj->M == 1) {
		for (i = 0; i < N; ++i) {
			obj->x[i] = obj->x[2 * N + i] - b[ip + iq];
		}
	}


	delta = 0.001;
	//mdisplay(P, 1, np);

	obj->ssq = ssq;
	//mdisplay(phi, 1, ip);
	ifault = flikam(phi, ip, theta, iq, obj->x, obj->x+N, N, &ssq, &fact, VW, VL, MRP1, VK, MR, delta);
	value = ssq*fact;

	//printf("sumlog ssq %g %g %d \n", sumlog,value,iter);

	free(phi);
	free(theta);
	free(VW);
	free(VL);
	free(VK);
	//free(bt);
	return value;
}

void as197(double *inp, int N, int optmethod, int p, int d, int q, double *phi, double *theta, double *wmean, double *var,double *resid,double *loglik,double *hess) {
	int i, pq, retval, length;
	double *b, *tf, *x;
	double maxstep;
	alik_object obj;
	//custom_function as154_min;

	x = (double*)malloc(sizeof(double)* (N - d));

	length = N;

	maxstep = 1.0;

	obj = alik_init(p, d, q, N);


	css(inp, N, optmethod, p, d, q, phi, theta, wmean, var,resid,loglik,hess);

	if (d > 0) {
		N = diff(inp, N, d, x); // No need to demean x
	}
	else {
		for (i = 0; i < N; ++i) {
			x[i] = inp[i];
		}
	}

	obj->N = N;
	obj->mean = *wmean;
	pq = obj->pq;
	b = (double*)malloc(sizeof(double)* pq);
	tf = (double*)malloc(sizeof(double)* pq);

	for (i = 0; i < p; ++i) {
		b[i] = phi[i];
	}
	for (i = 0; i < q; ++i) {
		b[p + i] = -theta[i];
	}

	if (obj->M == 1) {
		b[p + q] = obj->mean;
	}

	//mdisplay(b, 1, pq);

	for (i = 0; i < N; ++i) {
		obj->x[i] = obj->x[2 * N + i] = x[i];
	}
	for (i = N; i < 2 * N; ++i) {
		obj->x[i] = 0.0;
	}

	custom_function as197_min = { fas197, obj };
	retval = fminunc(&as197_min, NULL, pq, b, maxstep, optmethod, tf);

	for (i = 0; i < pq; ++i) {
		printf("%g ", tf[i]);
	}

	for (i = 0; i < p; ++i) {
		phi[i] = tf[i];
	}
	for (i = 0; i < q; ++i) {
		theta[i] = tf[p + i];
	}
	if (obj->M == 1) {
		*wmean = tf[p + q];
	}
	else {
		*wmean = 0.0;
	}
	for (i = 0; i < N - d; ++i) {
		resid[i] = obj->x[N + i];
	}
	*var = (obj->ssq) / (double)N;


	free(b);
	free(tf);
	free(x);
	free_alik(obj);
}

void free_alik_css(alik_css_object object) {
	free(object);
}

void free_alik(alik_object object) {
	free(object);
}

void free_alik_seas(alik_seas_object object) {
	free(object);
}

void free_alik_css_seas(alik_css_seas_object object) {
	free(object);
}
