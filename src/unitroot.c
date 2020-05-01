#include "unitroot.h"

void ur_df(double *y, int N,const char* type, int lags, const char *selectlags,double *statistic) {
    int lag,N1,i,j,N2,p,k,pN,tN,Nout;
    double *x,*z,*z_diff,*z_lags_1,*z_diff_lag,*res,*XX,*Y,*tt,*varcovar;
    reg_object fit;
    double alpha,ssr;
    double nT[6] = {25, 50, 100, 250, 500, 100000};
    double prob[8] ={0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99};
    double dftable[48] = {-4.38, -4.15, -4.04, -3.99, -3.98, -3.96,-3.95, -3.80, -3.73, -3.69, -3.68, -3.66,-3.60, -3.50, -3.45,
     -3.43, -3.42, -3.41,-3.24, -3.18, -3.15, -3.13, -3.13, -3.12,-1.14, -1.19, -1.22, -1.23, -1.24, -1.25,-0.80, -0.87, -0.90, -0.92,
      -0.93, -0.94,-0.50, -0.58, -0.62, -0.64, -0.65, -0.66,-0.15, -0.24, -0.28, -0.31, -0.32, -0.33};// 8X6
    double appxval[8] = {0,0,0,0,0,0,0,0};

    pN = 8;
    tN = 6;
    Nout = 50;

    if (lags < 0) {
        printf("Lags must be >= 0 \n");
        exit(-1);
    }

    alpha = 0.95;
    lag = lags;
    lags++;
    N1 = N - 1;
    N2 = N1 - lags + 1;

    z = (double*)malloc(sizeof(double)*N1);
    x = (double*)malloc(sizeof(double)*(N1-lags+1)*lags);
    res = (double*)malloc(sizeof(double)*N2);
    Y = (double*)malloc(sizeof(double)*N2);
    tt = (double*)malloc(sizeof(double)*N2);

    diff(y,N,1,z);// z = y(t) - y(t-1)

    
    for(i = 0; i < lags;++i) {
        for(j = 0;j < N2;++j) {
            x[i*N2+j] = z[lags+j-i-1];
        }
    } 

    mdisplay(x,lags,N1-lags+1);
    
    memcpy(Y,x,sizeof(double)*N2);

    z_diff = &x[0];// Length N2
    z_lags_1 = &y[lags-1];// Length N2 *(lags-1)

    for(i = 0; i < N2;++i) {
        tt[i] = (double) (i + lags);
    }
    

    if (lags > 1) {
        p = lags +1;
        XX = (double*)malloc(sizeof(double)*N2*p);
        varcovar = (double*)malloc(sizeof(double)*p*p);
        fit = reg_init(N2,p);
        memcpy(XX,x,sizeof(double)*N2);
        memcpy(XX+N2,tt,sizeof(double)*N2);
        memcpy(XX+2*N2,x+N2,sizeof(double)*N2*(lags-1));
        regress(fit,XX,Y,res,varcovar,alpha);
    } else {
        p = 2;
        XX = (double*)malloc(sizeof(double)*N2*p);
        varcovar = (double*)malloc(sizeof(double)*p*p);
        fit = reg_init(N2,p);
        memcpy(XX,x,sizeof(double)*N2);
        memcpy(XX+N2,tt,sizeof(double)*N2);
        regress(fit,XX,Y,res,varcovar,alpha);
    }

    *statistic = (fit->beta+1)->value / (fit->beta+1)->stdErr;

    for(i = 0; i < pN;++i) {

    }
    

    free(z);
    free(x);
    free(res);
    free(Y);
    free(tt);
    free(XX);
    free(varcovar);
    free_reg(fit);
}

void ur_kpss(double *y, int N,const char* type) {
    int i;
    double *tt;
    reg_object fit;

    tt = (double*)malloc(sizeof(double)*N);

    for(i = 1; i <= N;++i) {
        tt[i-1] = (double) i;
    }

    free(tt);

}
