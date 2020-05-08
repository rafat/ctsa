#include "unitroot.h"

void ur_df(double *y, int N,const char* alternative, int *klag, double *statistic,double *pval) {
    int lag,lags,N1,i,j,N2,p,k,pN,tN,Nout,iter;
    double *x,*z,*z_diff,*z_lags_1,*z_diff_lag,*res,*XX,*Y,*tt,*varcovar,*tablep;
    reg_object fit;
    double alpha,ssr,ylo,yhi,interp;
    double nT[6] = {25, 50, 100, 250, 500, 100000};
    double prob[8] ={0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99};
    double dftable[48] = {-4.38, -4.15, -4.04, -3.99, -3.98, -3.96,-3.95, -3.80, -3.73, -3.69, -3.68, -3.66,-3.60, -3.50, -3.45,
     -3.43, -3.42, -3.41,-3.24, -3.18, -3.15, -3.13, -3.13, -3.12,-1.14, -1.19, -1.22, -1.23, -1.24, -1.25,-0.80, -0.87, -0.90, -0.92,
      -0.93, -0.94,-0.50, -0.58, -0.62, -0.64, -0.65, -0.66,-0.15, -0.24, -0.28, -0.31, -0.32, -0.33};// 8X6

    pN = 8;
    tN = 6;
    Nout = 50;

    if (*klag < 0) {
        lags = (int) pow((double)N - 1.0, 1.0/3.0);
    } else {
        lags = *klag;
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
    tablep = (double*)malloc(sizeof(double)*pN);

    diff(y,N,1,z);// z = y(t) - y(t-1)

    //mdisplay(z,1,N1);
    
    for(i = 0; i < lags;++i) {
        for(j = 0;j < N2;++j) {
            x[i*N2+j] = z[lags+j-i-1];
        }
    } 

    //mdisplay(x,lags,N1-lags+1);
    
    memcpy(Y,x,sizeof(double)*N2);

    mdisplay(Y,1,N2);

    z_diff = &x[0];// Length N2
    z_lags_1 = &y[lags-1];// Length N2 *(lags-1)

    for(i = 0; i < N2;++i) {
        tt[i] = (double) (i + lags);
    }
    

    if (lags > 1) {
        //lags 5
        p = lags + 2;
        XX = (double*)malloc(sizeof(double)*N2*p);
        varcovar = (double*)malloc(sizeof(double)*p*p);
        fit = reg_init(N2,p);
        setLLSMethod(fit,"normal");
        memcpy(XX,y+lags-1,sizeof(double)*N2);
        memcpy(XX+N2,tt,sizeof(double)*N2);
        memcpy(XX+2*N2,x+N2,sizeof(double)*N2*(lags-1));
        mdisplay(XX,lags+1,N2);
        regress(fit,XX,Y,res,varcovar,alpha);
    } else {
        p = 3;
        XX = (double*)malloc(sizeof(double)*N2*p);
        varcovar = (double*)malloc(sizeof(double)*p*p);
        fit = reg_init(N2,p);
        memcpy(XX,y+lags-1,sizeof(double)*N2);
        memcpy(XX+N2,tt,sizeof(double)*N2);
        regress(fit,XX,Y,res,varcovar,alpha);
    }

    summary(fit);
    anova(fit);
    confint(fit);

    *statistic = (fit->beta+1)->value / (fit->beta+1)->stdErr;

    for(i = 0; i < pN;++i) {
        iter = i * tN;
        arrayminmax(dftable+iter,tN,&ylo,&yhi);
        tablep[i] = interpolate_linear(nT,dftable+iter,tN,ylo,yhi,(double)N1);
    }

    interp = interpolate_linear(tablep,prob,pN,prob[0],prob[pN-1],*statistic);

    //printf("interp %g \n",interp);

    if (!strcmp(alternative,"stationary")) {
        *pval = interp;
    } else if (!strcmp(alternative,"explosive")) {
        *pval = 1.0 - interp;
    } else {
        printf("alternative accepts only two values - stationary and explosive \n");
        exit(-1);
    }

    *klag = lags - 1;
    

    free(z);
    free(x);
    free(res);
    free(Y);
    free(tt);
    free(XX);
    free(varcovar);
    free(tablep);
    free_reg(fit);
}

void ur_kpss(double *y, int N,const char* type,int lshort,double *statistic,double *pval) {
    int i,p,l;
    double *tt,*res,*table,*varcovar,*csum;
    double alpha,eta,s2;
    reg_object fit;

    double tablep[4] = {0.01, 0.025, 0.05, 0.10};

    tt = (double*)malloc(sizeof(double)*N);
    res = (double*)malloc(sizeof(double)*N);
    table = (double*)malloc(sizeof(double)*4);
    csum = (double*)malloc(sizeof(double)*N);

    for(i = 1; i <= N;++i) {
        tt[i-1] = (double) i;
    }

    alpha = 0.95;

    if (!strcmp(type,"Trend")) {
            table[0] = 0.216; table[1] = 0.176; table[2] = 0.146; table[3] = 0.119;
            p = 1;
            varcovar = (double*)malloc(sizeof(double)*p*p);

            fit = reg_init(N,p);
            setIntercept(fit,0);
            regress(fit,tt,y,res,varcovar,alpha);

        } else if (!strcmp(type,"Level")) {
            table[0] = 0.739; table[1] = 0.574; table[2] = 0.463; table[3] = 0.347;
            p = 1;
            varcovar = (double*)malloc(sizeof(double)*p*p);

            fit = reg_init(N,p);
            regress(fit,NULL,y,res,varcovar,alpha);
    }

    cumsum(res,N,csum);

    eta = s2 = 0.0;

    for(i = 0; i < N;++i) {
        eta += (csum[i] * csum[i]);
        s2 += (res[i] * res[i]);
    }

    eta /= (double) (N*N);
    s2 /= (double) N;

    if (lshort == 1) {
        l = (int) 4.0 * pow(((double) N)/100.0,0.25);
    } else {
        l = (int) 12.0 * pow(((double) N)/100.0,0.25);
    }

    ppsum(res,N,l,&s2);

    *statistic = eta/s2;

    *pval = interpolate_linear(table,tablep,4,tablep[0],tablep[3],*statistic);

    free(table);
    free(tt);
    free(res);
    free(csum);
    free(varcovar);
    free_reg(fit);
}

void ur_pp(double *y, int N,const char* alternative,const char* type,int lshort) {
    int i,j,N1;
    double *z,*yt,*yt1,*tt;

    N1 = N - 1;

    z = (double*)malloc(sizeof(double)*2*N1);
    tt = (double*)malloc(sizeof(double)*N1);

    for(i = 0; i < 2;++i) {
        for(j = 0;j < N1;++j) {
            z[i*N1+j] = y[1+j-i];
        }
    } 

    yt = &z[0];
    yt1 = &z[N1];

    for(i = 1; i <= N1;++i) {
        tt[i-1] = (double)(i - N1);
        tt[i-1] /= 2.0;
    }

    free(z);
}
