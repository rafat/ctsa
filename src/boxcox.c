// SPDX-License-Identifier: BSD-3-Clause
#include "boxcox.h"

int checkConstant(double *x, int N) {
    int i, cc;
    double diff;
    cc = 0;

    /// heap buffer overflow fix, since accessing x[i+1] need to 
    /// stop loop at -1 index from end of array as dictated by length defined by N
    for(i = 1; i < N - 1;++i) {
        diff = x[i+1] - x[i];
        if (fabs(diff) > 0) {
            return cc;
        }
    }

    cc = 1;

    return cc;
}

static int checkNegative(double *x, int N) {
    int i,cn;
    cn = 1;

    for(i = 0; i < N;++i) {
        if (x[i] < 0) {
            return cn;
        }
    }

    cn = 0;

    return cn;
}

void boxcox_eval(double *x, int N, double lambda,double *bxcx) {
    double eps;
    int i;

    eps = macheps() / 744.44;

    if (fabs(lambda) < eps) {
        for(i = 0; i < N;++i) {
            bxcx[i] = log(x[i]);
        }
    } else {
        for(i = 0; i < N;++i) {
            bxcx[i] = expm1(lambda * log(x[i])) / lambda;
        }
    }

}

double inv_boxcox_eval(double *x,int N, double lambda,double *bxcx) {
    int i;

    if (lambda == 0) {
        for(i = 0; i < N;++i) {
            bxcx[i] = exp(x[i]);
        }
    } else {
        for(i = 0; i < N;++i) {
            bxcx[i] = exp(log1p(lambda * x[i])/lambda);
        }
    }
}

double boxcox_loglik(double lambda, void *params) {
    double loglik,variance,sum;
    int N,i;
    double *tdata;
    double *xpar = (double*) params;

    N = (int) xpar[0];
    tdata = (double*) malloc(sizeof(double) *N);

    sum = 0.0;
    for(i = 0; i < N;++i) {
        tdata[i] = log(xpar[i+1]);
        sum += tdata[i];
    }

    if (lambda == 0) {
        variance = var(tdata,N);
    } else {
        for(i = 0; i < N;++i) {
            tdata[i] = pow(xpar[i+1],lambda) / lambda;
        }
        variance = var(tdata,N);
    }

    loglik = (lambda-1)*sum - (xpar[0]/2.0)*log(variance);

    free(tdata);
    return -loglik;
}

double boxcox(double *x, int N,double *lambda,double *y) {
    double *params;
    double lmbd,a,b;
    int cc,cn;

    // Check if x has all constant values
    cc = checkConstant(x,N);

    if (cc == 1) {
        printf("Input is a constant vector. Exiting. \n");
        exit(-1);
    }

    // Check for negative values

    if (cn == 1) {
        printf("Input vector cannot have negative values. Exiting. \n");
        exit(-1);
    }

    params = (double*) malloc(sizeof(double)*(N+1));

    params[0] = (double) N;
    memcpy(params+1,x,sizeof(double)*N);

    custom_funcuni boxcox_mle = {boxcox_loglik,params};

    if (lambda == NULL) {
        a = -2;
        b = 2;
        lmbd = fminbnd(&boxcox_mle,a,b);
    } else {
        lmbd = *lambda;
    }

    boxcox_eval(x,N,lmbd,y);

    free(params);
    return lmbd;
}