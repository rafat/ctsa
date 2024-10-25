// SPDX-License-Identifier: BSD-3-Clause
/*
 * pred.c
 *
 *  Created on: Jul 14, 2014
 *      Author: Rafat Hussain
 */


#include "pred.h"

void predictarima(double *inp,int leninp,int p,int d, int q,double *phi, double *theta,
		double constant,int forlength,double *oup) {
	int lenzt,i,count,c2,K,L,lentheta,M,inilen,counter2,atlen2,ztc,atc,it;
	double *zt,*phi1,*theta1,*coeff,*coeff2,*ztl,*at,*delta,*sum1,*sum2;
	double temp,temp2;

	lenzt = leninp - d;

	zt = (double*) malloc(sizeof(double) * lenzt);
	phi1 = (double*) malloc(sizeof(double) * (p+1));
	theta1 = (double*) malloc(sizeof(double) * (q+1));
	delta = (double*) malloc(sizeof(double) * (d+1));
	at = (double*) malloc(sizeof(double) * (q+1+lenzt));
	ztl = (double*) malloc(sizeof(double) * (lenzt+forlength));


	if (d > 0) {
		lenzt = diff(inp,leninp,d,zt);
	} else {
		lenzt = leninp;
		for(i = 0; i < lenzt;++i) {
			zt[i] = inp[i];
		}
	}

	deld(d,delta);
	K = p+1;
	L = p + d + 1;
	lentheta = q + 1;

	sum1 = (double*) malloc(sizeof(double) * (d+1));
	sum2 = (double*) malloc(sizeof(double) * L);

	phi1[0] = theta1[0] = 1.0;

	for(i = 0; i < p;++i) {
		phi1[i+1] = - phi[i];
	}

	for(i = 0; i < q;++i) {
		theta1[i+1] = - theta[i];
	}

	for(i = 0; i < q+1+lenzt;++i) {
		at[i] = 0.0;
	}

	coeff = (double*) malloc(sizeof(double) * L);
	coeff2 = (double*) malloc(sizeof(double) * (L-1));

	//mdisplay(coeff2,1,L-1);

	//poly(phi1,delta,coeff,p+1,d+1);

	for(i = 0; i < L;++i) {
		coeff[i] = 0.0;
	}

	for(it = 0; it < K;++it) {
		for(i = 0; i < d+1;++i) {
			sum1[i] = delta[i] * phi1[it];
		}
		for(i = 0; i < L;++i) {
			sum2[i] = 0.0;
		}
		for(i = it; i < it+d+1;++i) {
			sum2[i] = sum1[i-it];
		}
		for(i = 0; i < L;++i) {
			coeff[i] += sum2[i];
		}

	}

	for(i = 0; i < L-1;++i) {
		coeff2[i] = -coeff[i+1];
	}


	inilen = 2;
	M = lenzt-inilen;

	counter2 = 0;
	atlen2 = lentheta;

	for(c2 = lenzt-inilen; c2 > 0; c2--) {

		for(i = 0; i < lenzt-c2+1;++i) {
			ztl[i] = zt[i];
		}

		for(count = 0; count < forlength;++count) {
			temp = constant;
			temp2 = 0.0;

			for(ztc = 0; ztc < L-1;++ztc) {
				temp=temp+coeff2[ztc]*ztl[lenzt-c2+count-ztc];
			}

			if (counter2 < M - 1) {
				for(atc = 0; atc < lentheta;++atc) {
					if ((atc > 0) && (count < lentheta-1) && (atlen2+count-atc <= atlen2-1)) {
						temp2=temp2+theta1[atc]*at[atlen2+count-atc];
					}
				}
			} else {
				for(atc = 0; atc < lentheta;++atc) {
					if ((atc > 0) && (count < lentheta-1) && (atlen2+count-atc <= atlen2-1)) {
						temp2=temp2+theta1[atc]*at[atlen2+count-atc];
					}
				}
			}
			ztl[lenzt-c2+1+count]=temp+temp2;
		}
		//mdisplay(ztl,1,lenzt+forlength);

		 if (c2 > 1) {
			 at[atlen2] = zt[counter2+inilen+1]-ztl[lenzt-c2+1];
			 atlen2++;
		 }
		 counter2++;
	}

	if (d == 0) {
		for(i = 0; i < forlength;++i) {
			oup[i] = ztl[lenzt+i];
		}
	} else {
		for (i = 0; i < forlength;++i) {
			oup[i] = ztl[lenzt+i];
			for(it = 1; it < d+1;++it) {
				if (leninp+i-it < leninp) {
					oup[i] -= (delta[it] * inp[leninp+i-it]);
				} else {
					oup[i] -= (delta[it] * oup[i-it]);
				}
			}
		}
	}

	free(zt);
	free(phi1);
	free(theta1);
	free(delta);
	free(coeff);
	free(coeff2);
	free(at);
	free(ztl);
	free(sum1);
	free(sum2);
}

