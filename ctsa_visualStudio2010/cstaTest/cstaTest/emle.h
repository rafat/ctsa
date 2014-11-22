/*
* emle.h
*
*  Created on: Jul 27, 2014
*      Author: Rafat Hussain
*/

#ifndef EMLE_H_
#define EMLE_H_

#include "initest.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct alik_css_set* alik_css_object;

alik_css_object alik_css_init(int p,int d,int q,int N);

struct alik_css_set{
	int p;// size of phi
	int d;// Number of times the series is to be differenced
	int q;//size of theta
	int r;// max(p,q+1)
	int pq;// p+q
	int length;// length of the original time series
	int N;// length of time series after differencing 
	int M;// M = 1 if mean needs to be calculated else 0
	double eps;
	double mean;
	double ssq;// Contains the sum of squares value 
	double loglik;
	double phi[100];
	double theta[100];
	double x[1];
};

typedef struct alik_set* alik_object;

alik_object alik_init(int p,int d,int q,int N);

struct alik_set{
	int p;// size of phi
	int d;// Number of times the series is to be differenced
	int q;//size of theta
	int r;// max(p,q+1)
	int pq;// p+q
	int length;// length of the original time series
	int N;// length of time series after differencing 
	int M;// M = 1 if mean needs to be calculated else 0
	double eps;
	double mean;
	double ssq;// Contains the sum of squares value 
	double loglik;
	double phi[100];
	double theta[100];
	double x[1];
};

typedef struct alik_css_seas_set* alik_css_seas_object;

alik_css_seas_object alik_css_seas_init(int p, int d, int q, int s, int P, int D, int Q, int N);

struct alik_css_seas_set{
	int p;// size of phi
	int d;// Number of times the series is to be differenced
	int q;//size of theta
	int s; // Frequency of Seasonal Components
	int P;// size of phi seasonal
	int D;// Number of times the seasonal series is to be differenced
	int Q;//size of theta seasonal
	int r;// max(p+s*P,q+s*Q+1)
	int pq;// p+q+s*P+s*Q
	int length;// length of the original time series
	int N;// length of time series after differencing 
	int M;// M = 1 if mean needs to be calculated else 0
	double eps;
	double mean;
	double ssq;// Contains the sum of squares value 
	double loglik;
	int offset;
	//double phi[100];
	//double theta[100];
	//double PHI[100];
	//double THETA[100];
	double x[1];
};

typedef struct alik_seas_set* alik_seas_object;

alik_seas_object alik_seas_init(int p, int d, int q, int s, int P, int D, int Q, int N);

struct alik_seas_set{
	int p;// size of phi
	int d;// Number of times the series is to be differenced
	int q;//size of theta
	int s; // Frequency of Seasonal Components
	int P;// size of phi seasonal
	int D;// Number of times the seasonal series is to be differenced
	int Q;//size of theta seasonal
	int r;// max(p+s*P,q+s*Q+1)
	int pq;// p+q+s*P+s*Q
	int length;// length of the original time series
	int N;// length of time series after differencing 
	int M;// M = 1 if mean needs to be calculated else 0
	double eps;
	double mean;
	double ssq;// Contains the sum of squares value 
	double loglik;
	int offset;
	//double phi[100];
	//double theta[100];
	//double PHI[100];
	//double THETA[100];
	double x[1];
};

int starma(int ip, int iq, double *phi, double *theta, double *A, double *P, double *V);

void karma(int ip, int iq, double *phi, double *theta, double *A, double *P, double*V, int N,
	double *W, double *resid, double *sumlog, double *ssq, int iupd, double delta, int *iter, int *nit);

int forkal(int ip, int iq, int id, double *phi, double*theta, double *delta, int N, double *W, double *resid, int il, double *Y, double *AMSE);

double fcss(double *b, int pq, void *params);

int css(double *inp, int N, int optmethod, int p, int d, int q, double *phi, double *theta, double *wmean, double *var,double *resid,double *loglik,double *hess);

double fcss_seas(double *b, int pq, void *params);

int css_seas(double *inp, int N, int optmethod, int p, int d, int q, int s, int P, int D, int Q,
	double *phi, double *theta, double *PHI, double *THETA, double *wmean, double *var,double *loglik,double *hess);

double fas154(double *x, int N, void *params);

int as154(double *x, int N, int optmethod, int p, int d, int q, double *phi, double *theta, double *wmean, double *var, double *resid, double *loglik, double *hess);

double fas154_seas(double *b, int pq, void *params);

int as154_seas(double *inp, int N, int optmethod, int p, int d, int q, int s, int P, int D, int Q,
	double *phi, double *theta, double *PHI, double *THETA, double *wmean, double *var,double *loglik,double *hess);

int flikam(double *P, int MP, double *Q, int MQ, double *W, double *E, int N, double *ssq, double *fact, double *VW, double *VL, int MRP1, double *VK, int MR, double TOLER);

double fas197(double *b, int pq, void *params);

void as197(double *inp, int N, int optmethod, int p, int d, int q, double *phi, double *theta, double *wmean, double *var,double *resid,double *loglik,double *hess);

void free_alik_css(alik_css_object object);

void free_alik(alik_object object);

void free_alik_seas(alik_seas_object object);

void free_alik_css_seas(alik_css_seas_object object);

#ifdef __cplusplus
}
#endif


#endif /* EMLE_H_ */