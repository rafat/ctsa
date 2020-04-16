#include "errors.h"

static void pe(double *predicted, double *actual, int N,double *err) {
    double temp;
    int i;


    for(i = 0; i < N;++i) {
        temp = (actual[i] - predicted[i]);

        err[i] = temp * 100 / actual[i];
    }

}

double me(double *predicted, double *actual, int N) {
    double err,temp;
    int i;

    temp = 0.0;

    for(i = 0; i < N;++i) {
        temp += (actual[i] - predicted[i]);
    }

    err = temp / N;

    return err;
}

double mse(double *predicted, double *actual, int N) {
    double err,temp;
    int i;

    temp = 0.0;

    for(i = 0; i < N;++i) {
        temp += ((actual[i] - predicted[i])*(actual[i] - predicted[i]));
    }

    err = temp / N;

    return err;
}

double rmse(double *predicted, double *actual, int N) {
    double err;

    err = mse(predicted,actual,N);

    err = sqrt(err);

    return err;
}

double mae(double *predicted, double *actual, int N) {
    double err,temp,t;
    int i;

    temp = 0.0;

    for(i = 0; i < N;++i) {
        t = fabs(actual[i] - predicted[i]);
        temp += t;
    }

    err = temp / N;

    return err;
}

double mape(double *predicted, double *actual, int N) {
    double err,temp;
    double *errvec;
    int i,j;

    errvec = (double*) malloc(sizeof(double)*N);

    pe(predicted,actual,N,errvec);
    j = 0;
    temp = 0.0;
    for(i = 0; i < N;++i) {
        if (errvec[i] == errvec[i]) {
            j++;
            temp += fabs(errvec[i]);
        }
    }

    err = temp/j;

    free(errvec);
    return err;
}

double mpe(double *predicted, double *actual, int N) {
    double err,temp;
    double *errvec;
    int i,j;

    errvec = (double*) malloc(sizeof(double)*N);

    pe(predicted,actual,N,errvec);
    j = 0;
    temp = 0.0;
    for(i = 0; i < N;++i) {
        if (errvec[i] == errvec[i]) {
            j++;
            temp += errvec[i];
        }
    }

    err = temp/j;

    free(errvec);
    return err;
}

double mase(double *predicted, double *actual, int N, double *tseries, int length) {
    double err,temp,den;
    double *dif;
    int i;

    // tseries : Only training values are included

    dif = (double*)malloc(sizeof(double)* (N-1));

    temp = 0.0;

    for(i = 1; i < length;++i) {
        temp += fabs(tseries[i] - tseries[i-1]);
    }

    den = temp / (length-1);

    temp = 0.0;

    for(i = 0; i < N;++i) {
        temp += fabs(actual[i] - predicted[i]);
    }

    err = temp /(N * den);


    return err;
}