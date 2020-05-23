#include "unitroot.h"

/*
 Modified C version of the R code by
 Adrian Trapletti and Kurt Hornik (2019). tseries: Time Series Analysis and Computational Finance. R package version 0.10-47. 
*/

void ur_df(double *y, int N,const char* alternative, int *klag, double *statistic,double *pval) {
    int lag,lags,N1,i,j,N2,p,k,pN,tN,Nout,iter;
    double *x,*z,*res,*XX,*Y,*tt,*varcovar,*tablep;
    reg_object fit;
    double alpha,ssr,interp;
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

    for(i = 0; i < N2;++i) {
        tt[i] = (double) (i + lags);
    }
    

    if (lags > 1) {
        //lags 5
        p = lags + 2;
        XX = (double*)malloc(sizeof(double)*N2*p);
        varcovar = (double*)malloc(sizeof(double)*p*p);
        fit = reg_init(N2,p);
        memcpy(XX,y+lags-1,sizeof(double)*N2);
        memcpy(XX+N2,tt,sizeof(double)*N2);
        memcpy(XX+2*N2,x+N2,sizeof(double)*N2*(lags-1));
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

    *statistic = (fit->beta+1)->value / (fit->beta+1)->stdErr;

    for(i = 0; i < pN;++i) {
        iter = i * tN;
        tablep[i] = interpolate_linear(nT,dftable+iter,tN,(double)N1);
    }

    interp = interpolate_linear(tablep,prob,pN,*statistic);

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

void ur_kpss(double *y, int N,const char* type,int lshort, int *klag, double *statistic,double *pval) {
    int i,p,l;
    double *tt,*res,*table,*varcovar,*csum;
    double alpha,eta,s2,ylo,yhi;
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
            p = 2;
            varcovar = (double*)malloc(sizeof(double)*p*p);

            fit = reg_init(N,p);
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

    arrayminmax(table,4,&ylo,&yhi);

    *pval = interpolate_linear(table,tablep,4,*statistic);

    *klag = l;

    free(table);
    free(tt);
    free(res);
    free(csum);
    free(varcovar);
    free_reg(fit);
}

void ur_pp(double *y, int N,const char* alternative,const char* type,int lshort, int *klag, double *statistic,double *pval) {
    int i,j,N1,N2,p,l,pN,tN,iter;
    double *z,*yt,*yt1,*tt,*XX,*varcovar,*res,*table,*tablep;
    double at,nssqr_res,trm1,trm2,trm3,trm4,dx,alpha,s2,tstat;
    double sum_yt1,sum_yt12,sum_yt1n,interp;
    long int N3;
    reg_object fit;

    double nT[6] = {25, 50, 100, 250, 500, 100000};
    double prob[8] ={0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99};
    double dftable_at[48] = {-4.38, -4.15, -4.04, -3.99, -3.98, -3.96,-3.95, -3.80, -3.73, -3.69, -3.68, -3.66,-3.60, -3.50, -3.45,
     -3.43, -3.42, -3.41,-3.24, -3.18, -3.15, -3.13, -3.13, -3.12,-1.14, -1.19, -1.22, -1.23, -1.24, -1.25,-0.80, -0.87, -0.90, -0.92,
      -0.93, -0.94,-0.50, -0.58, -0.62, -0.64, -0.65, -0.66,-0.15, -0.24, -0.28, -0.31, -0.32, -0.33};// 8X6

    double dftable[48] = {-22.5, -25.7, -27.4, -28.4, -28.9, -29.5, -19.9, -22.4, -23.6, -24.4, -24.8, -25.1, -17.9, -19.8, -20.7,
     -21.3, -21.5, -21.8, -15.6, -16.8, -17.5, -18.0, -18.1, -18.3, -3.66, -3.71, -3.74, -3.75, -3.76, -3.77, -2.51, -2.60, -2.62, 
     -2.64, -2.65, -2.66, -1.53, -1.66, -1.73, -1.78, -1.78, -1.79, -0.43, -0.65, -0.75, -0.82, -0.84, -0.87};// 8X6

    N1 = N - 1;
    N2 = N1 * N1;
    N3 = N2 * N1;

    pN = 8;
    tN = 6;

    z = (double*)malloc(sizeof(double)*2*N1);
    tt = (double*)malloc(sizeof(double)*N1);
    XX = (double*)malloc(sizeof(double)*2*N1);
    res = (double*)malloc(sizeof(double)*N1);
    tablep = (double*)malloc(sizeof(double)*pN);

    for(i = 0; i < 2;++i) {
        for(j = 0;j < N1;++j) {
            z[i*N1+j] = y[1+j-i];
        }
    } 

    yt = &z[0];
    yt1 = &z[N1];

    for(i = 1; i <= N1;++i) {
        tt[i-1] = (double)i - ((double)N1)/2.0;
    }

    p = 3;
    varcovar = (double*)malloc(sizeof(double)*p*p);

    memcpy(XX,tt,sizeof(double)*N1);
    memcpy(XX+N1,yt1,sizeof(double)*N1);

    at = 0.95;

    fit = reg_init(N1,p);

    //mdisplay(XX,2,N1);
    //mdisplay(yt,1,N1);

    regress(fit,XX,yt,res,varcovar,at);

    summary(fit);
    anova(fit);

    nssqr_res = 0.0;

    for(i =  0; i < N1;++i) {
        nssqr_res += (res[i]*res[i]);
    }

    nssqr_res /= ((double) N1);

    s2 = nssqr_res;

    if (lshort == 1) {
        l = (int) 4.0 * pow(((double) N1)/100.0,0.25);
    } else {
        l = (int) 12.0 * pow(((double) N1)/100.0,0.25);
    }

    ppsum(res,N1,l,&nssqr_res);

    sum_yt1 = sum_yt12 = sum_yt1n = 0.0;

    for(i = 0; i < N1;++i) {
        sum_yt1 += yt1[i];
        sum_yt12 += (yt1[i] * yt1[i]);
        sum_yt1n += (yt1[i] * (double) (i+1));
    }

    trm1 = N2 * (N2-1) *sum_yt12 / 12.0;
    trm2 = N1 * sum_yt1n * sum_yt1n;
    trm3 = N1 * (N1+1) * sum_yt1n * sum_yt1;
    trm4 = (N1 * (N1+1) * (2*N1+1) * sum_yt1 * sum_yt1) / 6.0;

    dx = trm1-trm2+trm3-trm4;

    printf("dx %g \n",dx);    

    if (!strcmp(type,"Z(alpha)")) {
        table = &dftable[0];
        alpha = (fit->beta+2)->value;
        *statistic = (double) N1 * (alpha - 1.0) - ((double)N3 / (24.0 * dx)) * ((double)N3 * (nssqr_res - s2));
    } else if (!strcmp(type,"Z(t_alpha)")) {
        table = &dftable_at[0];
        tstat = ((fit->beta+2)->value - 1.0) / (fit->beta+2)->stdErr;
        *statistic = sqrt(s2)/sqrt(nssqr_res)*tstat - ((double) N3) / (4*sqrt(3.0)*sqrt(dx)*sqrt(nssqr_res))*(nssqr_res - s2);
    } else {
        printf(" Type only admits two values : Z(alpha) and Z(t_alpha) \n");
        exit(-1);
    }

    printf("stat %g (fit->beta+2)->value %g \n",*statistic,(fit->beta+2)->value);

    for(i = 0; i < tN*pN;++i) {
        printf("%g ",table[i]);
    }

    printf("\n ss1 %g ss2 %g \n",s2,nssqr_res);

    *klag = l;

    
    for(i = 0; i < pN;++i) {
        iter = i * tN;
        tablep[i] = interpolate_linear(nT,table+iter,tN,(double)N1);
    }
    
    interp = interpolate_linear(tablep,prob,pN,*statistic);

    printf("interp %g \n",interp);

    if (!strcmp(alternative,"stationary")) {
        *pval = interp;
    } else if (!strcmp(alternative,"explosive")) {
        *pval = 1.0 - interp;
    } else {
        printf("alternative accepts only two values - stationary and explosive \n");
        exit(-1);
    }

    

    free_reg(fit);
    free(z);
    free(tt);
    free(XX);
    free(varcovar);
    free(res);
    free(tablep);
}
