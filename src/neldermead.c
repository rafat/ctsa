
#include "neldermead.h"

/*
 * neldermead.c
 *
 *  Copyright (c) 2014, Rafat Hussain
 *	License : BSD 3-Clause
 *	See COPYRIGHT for more details
 */
 
/*
 * C Routine nel_min is based on O'Neill's FORTRAN implementation of Nelder-Mead algorithm
 * Algorithm AS 47: Function Minimization Using a Simplex Procedure
	Author(s): R. O'Neill
	Source: Journal of the Royal Statistical Society. Series C (Applied Statistics), Vol. 20, No. 3
	(1971), pp. 338-345
	Published by: Blackwell

	Main Reference

	J. A. Nelder and R. Mead, A simplex method for function minimization. Compuer Journal,7,308-313

 */

int nel_min(custom_function *funcpt,double *xc,int N,double *dx,double fsval,int MAXITER,int *niter,
		double eps,double *xf) {
	int rcode,NN,i,j,ct,ctl;
	int ihi,ilo,L,siter;
	double reqmin,rcoeff,ecoeff,ccoeff,fn,del;
	double ylo,ynl,temp,ys,yss;
	double *P,*Y,*PB,*PS,*PSS,*PM,*PL,*xi,*xmin;
	double ysum,yavg,s,dval;

	rcode = 0;
	reqmin = eps;
	rcoeff = 1.0;
	ecoeff = 2.0;
	ccoeff = 0.5;
	del = 1.0;
	NN = N + 1;
	ct = N * N;
	*niter = 0;
	siter = MAXITER;
	s = 1.0;
	ctl = 0;
	dval = 1.0e-03;

	P = (double*) malloc(sizeof(double) * N * NN);
	Y = (double*) malloc(sizeof(double) * NN);
	PB = (double*) malloc(sizeof(double) * N);
	PS = (double*) malloc(sizeof(double) * N);
	PSS = (double*) malloc(sizeof(double) * N);
	PM = (double*) malloc(sizeof(double) * N);
	PL = (double*) malloc(sizeof(double) * N);
	xi = (double*) malloc(sizeof(double) * N);
	xmin = (double*) malloc(sizeof(double) * N);

	for (i = 0; i < N;++i) {
		xi[i] = xc[i];
	}

	while (ctl == 0 && *niter < siter) {
		// Create N+1 point Simplex
		for (i = 0; i < N;++i) {
			P[ct + i] = xi[i];
		}

		fn = FUNCPT_EVAL(funcpt,xi,N);
		Y[N] = fn;


		for (i = 0; i < N;++i) {
			xi[i] += dx[i] * del;
			ct = i * N;
			for (j = 0; j < N;++j) {
				P[ct + j] = xi[j];
			}
			fn = FUNCPT_EVAL(funcpt, xi, N);
			Y[i] = fn;
			xi[i] -= dx[i] * del;
		}

		ylo = Y[0];
		ilo = 0;

		for (i = 1; i < NN;++i) {

			if (Y[i] < ylo) {
				ylo = Y[i];
				ilo = i;
			}
		}

		ct = ilo * N;

		for (i = 0; i < N;++i) {
			PL[i] = P[ct + i];
		}


		while (fabs(s) > reqmin && *niter < siter) {
			*niter = *niter + 1;

			// Find highest and lowest function values.

			ynl = Y[0];
			ihi = 0;

			for (i = 1; i < NN;++i) {

				if (Y[i] > ynl) {
					ynl = Y[i];
					ihi = i;
				}

			}


			// Find PB [Centroid]
			ct = N;
			for(i = 0; i < N;i++) {
				temp = 0.0;
				for ( j = 0; j < NN;++j) {
					temp += P[j*ct+i];
				}
				temp -= P[ihi*ct+i];
				PB[i] = temp / (double) N;
			}

			// Find reflection
			ct = ihi * N;
			for (i = 0; i < N;++i) {
				PS[i] = (1.0 + rcoeff) * PB[i] - rcoeff * P[ct+i];
			}
			ys = FUNCPT_EVAL(funcpt, PS, N);

			if (ys < ylo) {
				for (i = 0; i < N;++i) {
					PSS[i] = ecoeff * PS[i] + (1.0 - ecoeff) * PB[i];
				}
				yss = FUNCPT_EVAL(funcpt, PSS, N);

				if (yss < ys) {
					for(i = 0; i < N;++i) {
						P[ct + i] = PSS[i];
					}
					Y[ihi] = yss;
				} else {
					for(i = 0; i < N;++i) {
						P[ct + i] = PS[i];
					}
					Y[ihi] = ys;
				}


			} else {
				L = 0;
				for (i = 0; i < NN;++i) {
					if (Y[i] > ys) {
						L++;
					}
				}
	/*
				if (L == 1) {
					for ( i = 0; i < N;++i) {
						P[ct + i] = PS[i];
					}
					Y[ihi] = ys;
				}*/
				if (L == 0) {
					for (i = 0; i < N;++i) {
						PSS[i] = ccoeff * P[ct + i] + (1.0 - ccoeff) * PB[i];
					}
					yss = FUNCPT_EVAL(funcpt, PSS, N);

					if (yss <= Y[ihi]) {
						for(i = 0; i < N;++i) {
							P[ct + i] = PSS[i];
						}
						Y[ihi] = yss;

					} else {

						//ctl = ilo * N;
						for(i = 0; i < NN;++i) {
							ct = i * N;
							for (j = 0;j < N;++j) {
								P[ct + j] = 0.5 * (P[ct + j] + PL[j]);
								PM[j] = P[ct + j];
							}
							Y[i] = FUNCPT_EVAL(funcpt, PM, N);

						}
						ylo = Y[0];
						ilo = 0;

						for (i = 1; i < NN;++i) {

							if (Y[i] < ylo) {
								ylo = Y[i];
								ilo = i;
							}

						}
						ct = ilo * N;

						for (i = 0; i < N;++i) {
							PL[i] = P[ct + i];
						}
					}
				}
				if (L == 1) {
					for (i = 0; i < N;++i) {
						PSS[i] = ccoeff * P[ct + i] + (1.0 - ccoeff) * PB[i];
					}
					yss = FUNCPT_EVAL(funcpt, PSS, N);

					if (yss <= ys) {
						ct = ihi * N;
						for(i = 0; i < N;++i) {
							P[ct + i] = PSS[i];
						}
						Y[ihi] = yss;

					} else {
						ct = ihi * N;
						for(i = 0; i < N;++i) {
							P[ct + i] = PS[i];
						}
						Y[ihi] = ys;

					}

				}


				if (L > 1) {
					ct = ihi * N;
					for(i = 0; i < N;++i) {
						P[ct + i] = PS[i];
					}
					Y[ihi] = ys;

				}
			}

			if (Y[ihi] < ylo) {
				ylo = Y[ihi];
				ilo = ihi;
			}

			//To the end of the loop
			ysum = 0.0;
			for (j=0 ;j <= N;j++) {
				ysum += Y[j];
			}
			yavg = ysum/(N+1);
			s = 0.0;
			for (j=0;j<=N;j++) {
				s += pow((Y[j]-yavg),2.0);
			}
			s = sqrt(s);//Remove this statement if you want convergence in fewer iterations
			if (fabs(s) < reqmin) {
				ctl = 1;
			}

		}
		ct = ilo * N;
		for(i = 0; i < N;++i) {
			xmin[i] = P[ct + i];
		}
		ynl = Y[ilo];

		for (i = 0; i < N;++i) {
			del = dx[i] * dval;
			xmin[i] += del;
			fn = FUNCPT_EVAL(funcpt, xmin, N);
			if (fn < ynl) {
				ctl = 0;
				break;
			} else {
				ctl++;
			}
			xmin[i] = xmin[i] - 2 * del;
			fn = FUNCPT_EVAL(funcpt, xmin, N);
			if (fn < ynl) {
				ctl = 0;
				break;
			} else {
				ctl++;
			}
		}

		if (ctl == 0) {
			for (i = 0; i < N;++i) {
				xi[i] = xmin[i];
			}
			del = dval;
			s = 1.0;
		}

	}

	ct = N * N;
	for (i = 0; i < N;++i) {
		xf[i] = P[ct + i];
	}

	if (rcode == 0 && *niter >= siter) {
		rcode = 4;
	} else if (*niter < siter) {
		rcode = 1;
	}


	free(Y);
	free(P);
	free(PB);
	free(PS);
	free(PSS);
	free(PM);
	free(PL);
	free(xi);
	free(xmin);
	return rcode;
}
