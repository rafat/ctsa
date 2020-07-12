/*
* arma.h
*
*  Created on: Sep 11, 2014
*      Author: Rafat Hussain
*/

#ifndef CTSA_H_
#define CTSA_H_

#include "emle.h"
#include "autoutils.h"

#ifdef __cplusplus
extern "C" {
#endif
	/*
	ARIMA vs SARIMA
	Typically you should use SARIMA if you are calculating seasonal models and ARIMA for non-seasonal models. However, ARIMA wil not work if the number of parameters p and q exceed 100.
	For these extreme-case non-seasonal models, use SARIMA and set seasonal paramters to zero. SARIMA and ARIMA work identically for non-seasonal models otherwise.
	*/

typedef struct arima_set* arima_object;

arima_object arima_init(int p, int d, int q, int N);

struct arima_set{
	int N;// length of time series
	int Nused;//length of time series after differencing, Nused = N - d
	int method;
	int optmethod;
	int p;// size of phi
	int d;// Number of times the series is to be differenced
	int q;//size of theta
	int M; // M = 0 if mean is 0.0 else M = 1
	int ncoeff;// Total Number of Coefficients to be estimated
	int cssml;// Uses CSS before MLE if 1 else uses MLE only if cssml is 0
	double *phi;
	double *theta;
	double *vcov;// Variance-Covariance Matrix Of length lvcov
	int lvcov; //length of VCOV
	double *res;
	double mean;
	double var;
	double loglik;
	double aic;
	int retval;
	double params[0];
};

typedef struct sarima_set* sarima_object;

sarima_object sarima_init(int p, int d, int q,int s,int P,int D,int Q, int N);

struct sarima_set{
	int N;// length of time series
	int Nused;//length of time series after differencing, Nused = N - d - s*D
	int method;
	int optmethod;
	int p;// size of phi
	int d;// Number of times the series is to be differenced
	int q;//size of theta
	int s;// Seasonality/Period
	int P;//Size of seasonal phi
	int D;// The number of times the seasonal series is to be differenced
	int Q;//size of Seasonal Theta
	int M; // M = 0 if mean is 0.0 else M = 1
	int ncoeff;// Total Number of Coefficients to be estimated
	int cssml;// Uses CSS before MLE if 1 else uses MLE only if cssml is 0
	double *phi;
	double *theta;
	double *PHI;
	double *THETA;
	double *vcov;// Variance-Covariance Matrix Of length lvcov
	int lvcov; //length of VCOV
	double *res;
	double mean;
	double var;
	double loglik;
	double aic;
	int retval;
	double params[0];
};

typedef struct ar_set* ar_object;

ar_object ar_init(int method, int N);

struct ar_set{
	int N;// length of time series
	int method;
	int optmethod;// Valid only for MLE estimation
	int p;// size of phi
	int order; // order = p
	int ordermax; // Set Maximum order to be fit
	double *phi;
	double *res;
	double mean;
	double var;
	double aic;
	int retval;
	double params[0];
};

typedef struct sarimax_set* sarimax_object;

sarimax_object sarimax_init(int p, int d, int q,int P, int D, int Q,int s, int r,int imean, int N);

struct sarimax_set{
	int N;// length of time series
	int Nused;//length of time series after differencing, Nused = N - d
	int method;
	int optmethod;
	int p;// size of phi
	int d;// Number of times the series is to be differenced
	int q;//size of theta
	int s;// Seasonality/Period
	int P;//Size of seasonal phi
	int D;// The number of times the seasonal series is to be differenced
	int Q;//size of Seasonal Theta
	int r;// Number of exogenous variables
	int M; // M = 0 if mean is 0.0 else M = 1
	int ncoeff;// Total Number of Coefficients to be estimated
	double *phi;
	double *theta;
	double *PHI;
	double *THETA;
	double *exog;
	double *vcov;// Variance-Covariance Matrix Of length lvcov
	int lvcov; //length of VCOV
	double *res;
	double mean;
	double var;
	double loglik;
	double aic;
	int retval;
	int start;
	int imean;
	double params[0];
};

typedef struct auto_arima_set* auto_arima_object;

auto_arima_object auto_arima_init(int *pdqmax,int *PDQmax,int s, int r, int N);

struct auto_arima_set{
	int N;// length of time series
	int Nused;//length of time series after differencing, Nused = N - d
	int method;
	int optmethod;
	int pmax;// Maximum size of phi
	int dmax;// Maximum Number of times the series is to be differenced
	int qmax;// Maximum size of theta
	int Pmax;//Maximum Size of seasonal phi
	int Dmax;// Maximum number of times the seasonal series is to be differenced
	int Qmax;//Maximum size of Seasonal Theta
	int p;// size of phi
	int d;// Number of times the series is to be differenced
	int q;//size of theta
	int s;// Seasonality/Period
	int P;//Size of seasonal phi
	int D;// The number of times the seasonal series is to be differenced
	int Q;//size of Seasonal Theta
	int r;// Number of exogenous variables
	int M; // M = 0 if mean is 0.0 else M = 1
	int ncoeff;// Total Number of Coefficients to be estimated
	double *phi;
	double *theta;
	double *PHI;
	double *THETA;
	double *exog;
	double *vcov;// Variance-Covariance Matrix Of length lvcov
	int lvcov; //length of VCOV
	double *res;
	double mean;
	double var;
	double loglik;
	double ic;
	int retval;
	int start;
	int imean;
	int idrift;
	int stationary;
	int seasonal;
	int Order_max;
	int p_start;
	int q_start;
	int P_start;
	int Q_start;
	char information_criteria[10];
	int stepwise;
	int num_models;
	int approximation;
	char test[10];
	char type[10];
	char seas[10];
	double alpha_test;
	double alpha_seas;
	double lambda;
	double params[0];
};


typedef struct sarimax_wrapper_set* sarimax_wrapper_object;

sarimax_wrapper_object sarimax_wrapper(sarimax_wrapper_object model,double *y, int N,int *order, int *seasonal, double *xreg, int r, int drift,int mean,
	double *lambda, int biasadj,int method);

struct sarimax_wrapper_set{
	sarimax_object sarimax;
	int idrift;
	double aic;
	double bic;
	double aicc;
	double lambda;
	double sigma2;
};

typedef struct myarima_set* myarima_object;

myarima_object myarima(double *x, int N, int *order, int *seasonal, int constant, const char* ic, int trace, int approx,
	double offset, double *xreg, int r, int *method);

struct myarima_set{
	sarimax_object sarimax;
	int idrift;
	double ic;
	double aic;
	double bic;
	double aicc;
	double sigma2;
};

typedef struct aa_ret_set* aa_ret_object;

aa_ret_object auto_arima1(double *y, int N, int *ordermax, int *seasonalmax,int *maxcoeff, int s,int *DD, int *dd, int *start, int *stationary, int *seasonal, 
	const char *ic, int *stepwise, int *nmodels,int *approximation,int *method,double *xreg, int r, const char *test,const char *type, double *test_alpha, 
	const char *seas, double *seas_alpha, int *allowdrift, int *allowmean, double *lambda);

struct aa_ret_set{
	sarimax_wrapper_object Arima;
	myarima_object myarima;
	int otype;
};

myarima_object search_arima(double *x, int N,int d, int D, int p_max, int q_max, int P_max, int Q_max, int Order_max, int stationary,int s, const char *ic,
	int approximation, double *xreg, int r, double offset,int allowdrift, int allowmean, int method);

void arima_exec(arima_object obj, double *x);

void sarimax_exec(sarimax_object obj, double *inp,double *xreg) ;

void sarima_exec(sarima_object obj, double *x);

void auto_arima_exec(auto_arima_object obj, double *inp,double *xreg);

void ar_exec(ar_object obj, double *inp);

void arima_predict(arima_object obj, double *inp, int L, double *xpred, double *amse);

void sarima_predict(sarima_object obj, double *inp, int L, double *xpred, double *amse);

void sarimax_predict(sarimax_object obj, double *inp, double *xreg, int L,double *newxreg, double *xpred, double *amse);

void sarimax_wrapper_predict(sarimax_wrapper_object obj, double *inp, double *xreg, int L,double *newxreg, double *xpred, double *amse);

void ar_predict(ar_object obj, double *inp, int L, double *xpred, double *amse);

void ar(double *inp, int N, int p, int method, double *phi,double *var);

void arima_setMethod(arima_object obj, int value);

void sarima_setMethod(sarima_object obj, int value);

void sarimax_setMethod(sarimax_object obj, int value);

void arima_setCSSML(arima_object obj, int cssml);

void sarima_setCSSML(sarima_object obj, int cssml);

void arima_setOptMethod(arima_object obj, int value);

void sarima_setOptMethod(sarima_object obj, int value);

void sarimax_setOptMethod(sarimax_object obj, int value);

void arima_vcov(arima_object obj, double *vcov);

void sarima_vcov(sarima_object obj, double *vcov);

void sarimax_vcov(sarimax_object obj, double *vcov);

void arima_summary(arima_object obj);

void sarimax_summary(sarimax_object obj);

void sarimax_wrapper_summary(sarimax_wrapper_object obj);

void aa_ret_summary(aa_ret_object obj);

void sarima_summary(sarima_object obj);

void myarima_summary(myarima_object obj);

int ar_estimate(double *x, int N, int method);

void ar_summary(ar_object obj);

void model_estimate(double *x, int N, int d, int pmax, int h);

void pacf(double* vec, int N, double* par, int M);

void pacf_opt(double* vec, int N, int method, double* par, int M);

void acvf(double* vec, int N, double* par, int M);

void acvf_opt(double* vec, int N, int method, double* par, int M);

void acvf2acf(double *acf, int M);

void arima_free(arima_object object);

void sarimax_free(sarimax_object object);

void sarimax_wrapper_free(sarimax_wrapper_object object);

void myarima_free(myarima_object object);

void sarima_free(sarima_object object);

void aa_ret_free(aa_ret_object object);

void ar_free(ar_object object);

// Yule-Walker, Burg and Hannan Rissanen Algorithms for Initial Parameter Estimation

void yw(double *x, int N, int p, double *phi, double *var);

void burg(double *x, int N, int p, double *phi, double *var);

void hr(double *x, int N, int p, int q, double *phi, double *theta, double *var);


#ifdef __cplusplus
}
#endif


#endif /* CTSA_H_ */
