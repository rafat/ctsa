/*
* arma.h
*
*  Created on: Sep 11, 2014
*      Author: Rafat Hussain
*/

#ifndef CTSA_H_
#define CTSA_H_

#ifdef __cplusplus
extern "C" {
#endif

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
	double sigma2;
	double aic;
	double bic;
	double aicc;
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

void sarimax_exec(sarimax_object obj, double *inp,double *xreg) ;

void arima_exec(arima_object obj, double *x);

void sarima_exec(sarima_object obj, double *x);

void auto_arima_exec(auto_arima_object obj, double *inp,double *xreg);

void ar_exec(ar_object obj, double *inp);

void arima_predict(arima_object obj, double *inp, int L, double *xpred, double *amse);

void sarima_predict(sarima_object obj, double *inp, int L, double *xpred, double *amse);

void sarimax_predict(sarimax_object obj, double *inp, double *xreg, int L,double *newxreg, double *xpred, double *amse);

void auto_arima_predict(auto_arima_object obj, double *inp, double *xreg, int L,double *newxreg, double *xpred, double *amse);

void ar_predict(ar_object obj, double *inp, int L, double *xpred, double *amse);

void ar(double *inp, int N, int p, int method, double *phi,double *var);

void arima_setMethod(arima_object obj, int value);

void sarima_setMethod(sarima_object obj, int value);

void auto_arima_setMethod(auto_arima_object obj, int value);

void sarimax_setMethod(sarimax_object obj, int value);

void arima_setOptMethod(arima_object obj, int value);

void sarima_setOptMethod(sarima_object obj, int value);

void sarimax_setOptMethod(sarimax_object obj, int value);

void auto_arima_setOptMethod(auto_arima_object obj, int value);

void arima_vcov(arima_object obj, double *vcov);

void sarima_vcov(sarima_object obj, double *vcov);

void sarimax_vcov(sarimax_object obj, double *vcov);

void auto_arima_setApproximation(auto_arima_object obj, int approximation);

void auto_arima_setStepwise(auto_arima_object obj, int stepwise);

void auto_arima_setStationary(auto_arima_object obj, int stationary);

void auto_arima_setSeasonal(auto_arima_object obj, int seasonal);

void auto_arima_setStationarityParameters(auto_arima_object obj,const char *test, double alpha, const char *type);

void auto_arima_setSeasonalParameters(auto_arima_object obj,const char *test, double alpha);

void arima_summary(arima_object obj);

void sarima_summary(sarima_object obj);

void sarimax_summary(sarimax_object obj);

void auto_arima_summary(auto_arima_object obj);

int ar_estimate(double *x, int N, int method);

void ar_summary(ar_object obj);

void model_estimate(double *x, int N, int d, int pmax, int h);

void pacf(double* vec, int N, double* par, int M);

void pacf_opt(double* vec, int N, int method, double* par, int M);

void acvf(double* vec, int N, double* par, int M);

void acvf_opt(double* vec, int N, int method, double* par, int M);

void acvf2acf(double *acf, int M);

void arima_free(arima_object object);

void sarima_free(sarima_object object);

void sarimax_free(sarimax_object object);

void auto_arima_free(auto_arima_object object);

void ar_free(ar_object object);

// Yule-Walker, Burg and Hannan Rissanen Algorithms for Initial Parameter Estimation

void yw(double *x, int N, int p, double *phi, double *var);

void burg(double *x, int N, int p, double *phi, double *var);

void hr(double *x, int N, int p, int q, double *phi, double *theta, double *var);


#ifdef __cplusplus
}
#endif


#endif /* CTSA_H_ */
