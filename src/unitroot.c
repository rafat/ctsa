// SPDX-License-Identifier: BSD-3-Clause
#include "unitroot.h"

/*
 Modified C version of the R code by
 Adrian Trapletti and Kurt Hornik (2019). tseries: Time Series Analysis and Computational Finance. R package version 0.10-47. 
*/

/* 
 ur_pp2
 C version of the R code by
 Bernhard Pfaff [aut, cre], Eric Zivot [ctb], Matthieu Stigler [ctb]. urca: Unit Root and Cointegration Tests for Time Series Data
 License Information : https://cran.r-project.org/web/packages/urca/urca.pdf
*/

void ur_df(double *y, int N,const char* alternative, int *klag, double *statistic,double *pval) {
    /*
    y - Time Series data of length N
    alternative - "stationary" or "explosive"

    klag - Length of the lag. Uses pow((double)N - 1.0, 1.0/3.0) if klag == NULL
    Outputs

    statistic - Test statistic
    pval - p-value used to determine null hypothesis

    */
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

    lags = (klag == NULL) ? (int) pow((double)N - 1.0, 1.0/3.0) : *klag;

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

    if (klag != NULL) *klag = lags - 1;
    

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

void ur_df2(double *y, int N,const char* type, int *lags,const char *selectlags,double *cval,int *cvrows, int *cvcols, double *cprobs, double *teststat,int *ltstat) {
    int lags_, lag, N1, N2,i,j,iter,p,ltmp,p1,p2,p3,rowselec,row1;
    double *z, *x,*z_diff,*z_lag_1,*tt,*critRes,*z_diff_lag,*XX,*varcovar,*res,*XX2,*XX3;
    reg_object fit, phi1_fit,phi2_fit,phi3_fit;
    double alpha = 0.95;
    double ctemp,tau,scale,sos,dfs,phi1,phi2,phi3;
    double cval_tau1[18] = {-2.66, -1.95, -1.60,-2.62, -1.95, -1.61,-2.60, -1.95, -1.61,-2.58, -1.95, -1.62,
        -2.58, -1.95, -1.62,-2.58, -1.95, -1.62};
    double cval_tau2[18] = {-3.75, -3.00, -2.63,-3.58, -2.93, -2.60,-3.51, -2.89, -2.58,-3.46, -2.88, -2.57,
        -3.44, -2.87, -2.57,-3.43, -2.86, -2.57};
    double cval_tau3[18] = {-4.38, -3.60, -3.24,-4.15, -3.50, -3.18,-4.04, -3.45, -3.15,-3.99, -3.43, -3.13,
        -3.98, -3.42, -3.13,-3.96, -3.41, -3.12};
    double cval_phi1[18] = {7.88, 5.18, 4.12,7.06, 4.86, 3.94,6.70, 4.71, 3.86,6.52, 4.63, 3.81,6.47, 4.61,
         3.79,6.43, 4.59, 3.78};
    double cval_phi2[18] = {8.21, 5.68, 4.67,7.02, 5.13, 4.31,6.50, 4.88, 4.16,6.22, 4.75, 4.07,6.15, 4.71,
         4.05,6.09, 4.68, 4.03};
    double cval_phi3[18] = {10.61, 7.24, 5.91,9.31, 6.73, 5.61,8.73, 6.49, 5.47,8.43, 6.49, 5.47,8.34, 6.30,
         5.36,8.27, 6.25, 5.34};

    lags_ = (lags == NULL) ? 1 : *lags;
    lag = lags_;

    if (lag < 0) {
        printf("lags should be > 0. \n");
        exit(-1);
    }

    lags_+=1;

    N1 = N - 1;
    N2 = N1 - lags_ + 1;

    z = (double*)malloc(sizeof(double)*N1);
    x = (double*)malloc(sizeof(double)*N2*lags_);
    

    diff(y,N,1,z);// z = y(t) - y(t-1)
    
    for(i = 0; i < lags_;++i) {
        for(j = 0;j < N2;++j) {
            x[i*N2+j] = z[lags_+j-i-1];
        }
    } 


    z_diff = &x[0];//length N2
    z_lag_1 = &y[lags_-1];//length N2


    tt = (double*)malloc(sizeof(double)*N2);
    res = (double*)malloc(sizeof(double)*N2);

    for(i = lags_ - 1; i < N1;++i) {
        tt[i - lags_ +1] = i+1;
    }

    if (lags_ > 1) {
        if (strcmp(selectlags,"fixed")) {
            critRes = (double*) malloc(sizeof(double)*lags_);
            XX = (double*)malloc(sizeof(double)*N2*(lags_+1));
            for(i = 0; i < lags_;++i) {
                critRes[i] = NAN;
            }

            for(i = 1; i < lags_;++i) {
                iter = i * N2;

                memcpy(XX,z_lag_1,sizeof(double)*N2);

                if (!strcmp(type,"none")) {
                    p = 1 + i;
                    varcovar = (double*)malloc(sizeof(double)*p*p);
                    fit = reg_init(N2,p);
                    setIntercept(fit,0);
                    memcpy(XX+N2,x+N2,sizeof(double)*iter);
                    regress(fit,XX,z_diff,res,varcovar,alpha);
                    
                } else if (!strcmp(type,"drift")) {
                    p = 2 + i;
                    varcovar = (double*)malloc(sizeof(double)*p*p);
                    fit = reg_init(N2,p);
                    setIntercept(fit,1);
                    memcpy(XX+N2,x+N2,sizeof(double)*iter);
                    regress(fit,XX,z_diff,res,varcovar,alpha);

                } else if (!strcmp(type,"trend")) {
                    p = 3 + i;
                    varcovar = (double*)malloc(sizeof(double)*p*p);
                    fit = reg_init(N2,p);
                    setIntercept(fit,1);
                    memcpy(XX+N2,tt,sizeof(double)*N2);
                    memcpy(XX+2*N2,x+N2,sizeof(double)*iter);
                    regress(fit,XX,z_diff,res,varcovar,alpha);

                } else {
                    printf("type only accepts one of three values - none, drift and trend \n");
                    exit(-1);
                }

                if (!strcmp(selectlags,"aic")) {
                    critRes[i] = fit->aic;
                } else if (!strcmp(selectlags,"bic")) {
                    critRes[i] = fit->bic;
                }

                free_reg(fit);
            }

            ctemp = DBL_MAX;
            ltmp = 0;

            for(i = 1; i < lags_;++i) {
                if (critRes[i] < ctemp) {
                    ctemp = critRes[i];
                    ltmp = i;
                }
            }

            lags_ = ltmp+1;


            free(critRes);
            free(XX);
            free(varcovar);
        }

        XX = (double*)malloc(sizeof(double)*N2*(lags_+1));

        memcpy(XX,z_lag_1,sizeof(double)*N2);
        iter = (lags_-1) * N2;

        if (!strcmp(type,"none")) {
            p = lags_;
            varcovar = (double*)malloc(sizeof(double)*p*p);
            fit = reg_init(N2,p);
            setIntercept(fit,0);
            memcpy(XX+N2,x+N2,sizeof(double)*iter);
            regress(fit,XX,z_diff,res,varcovar,alpha);
            tau = (fit->beta + 0)->value / (fit->beta + 0)->stdErr;
            teststat[0] = tau;
            *ltstat = 1;
        } else if (!strcmp(type,"drift")) {
            p = 1 + lags_;
            varcovar = (double*)malloc(sizeof(double)*p*p);
            fit = reg_init(N2,p);
            setIntercept(fit,1);
            memcpy(XX+N2,x+N2,sizeof(double)*iter);
            regress(fit,XX,z_diff,res,varcovar,alpha);
            tau = (fit->beta + 1)->value / (fit->beta + 1)->stdErr;

            scale = fit->RSS/(double)fit->df_RSS;

            XX2 = (double*)malloc(sizeof(double)*N2*(lags_-1));
            memcpy(XX2,x+N2,sizeof(double)*N2*(lags_-1));
            p1 = lags_ - 1;
            phi1_fit = reg_init(N2,p1);
            setIntercept(phi1_fit,0);
            regress(phi1_fit,XX2,z_diff,res,varcovar,alpha);

            sos = phi1_fit->RSS - fit->RSS;
            dfs = (double) (phi1_fit->df_RSS - fit->df_RSS);

            phi1 = sos/dfs/scale;


            teststat[0] = tau;
            teststat[1] = phi1;
            *ltstat = 2;

            free(XX2);
            free(phi1_fit);
        } else if (!strcmp(type,"trend")) {
            p = 2 + lags_;
            varcovar = (double*)malloc(sizeof(double)*p*p);
            fit = reg_init(N2,p);
            setIntercept(fit,1);
            memcpy(XX+N2,tt,sizeof(double)*N2);
            memcpy(XX+2*N2,x+N2,sizeof(double)*iter);
            regress(fit,XX,z_diff,res,varcovar,alpha);
            tau = (fit->beta + 1)->value / (fit->beta + 1)->stdErr;

            scale = fit->RSS/(double)fit->df_RSS;

            XX2 = (double*)malloc(sizeof(double)*N2*(lags_-1));
            memcpy(XX2,x+N2,sizeof(double)*N2*(lags_-1));
            p2 = lags_ - 1;
            phi2_fit = reg_init(N2,p2);
            setIntercept(phi2_fit,0);
            regress(phi2_fit,XX2,z_diff,res,varcovar,alpha);

            sos = phi2_fit->RSS - fit->RSS;
            dfs = (double) (phi2_fit->df_RSS - fit->df_RSS);

            phi2 = sos/dfs/scale;


            free(XX2);
            free(phi2_fit);

            XX3 = (double*)malloc(sizeof(double)*N2*(lags_-1));
            memcpy(XX3,x+N2,sizeof(double)*N2*(lags_-1));
            p3 = lags_;
            phi3_fit = reg_init(N2,p3);
            setIntercept(phi3_fit,1);
            regress(phi3_fit,XX3,z_diff,res,varcovar,alpha);

            sos = phi3_fit->RSS - fit->RSS;
            dfs = (double) (phi3_fit->df_RSS - fit->df_RSS);

            phi3 = sos/dfs/scale;

            teststat[0] = tau;
            teststat[1] = phi2;
            teststat[2] = phi3;
            *ltstat = 3;

            free(XX3);
            free(phi3_fit);
        } else {
            printf("type only accepts one of three values - none, drift and trend \n");
            exit(-1);
        }

        free(XX);
        free(varcovar);

    } else {
        XX = (double*)malloc(sizeof(double)*2*N2);
        if (!strcmp(type,"none")) {
            p = 1;
            varcovar = (double*)malloc(sizeof(double)*p*p);
            fit = reg_init(N2,p);
            setIntercept(fit,0);
            memcpy(XX,z_lag_1,sizeof(double)*N2);
            regress(fit,XX,z_diff,res,varcovar,alpha);
            tau = (fit->beta + 0)->value / (fit->beta + 0)->stdErr;
            teststat[0] = tau;
            *ltstat = 1;
        } else if (!strcmp(type,"drift")) {
            p = 2;
            varcovar = (double*)malloc(sizeof(double)*p*p);
            fit = reg_init(N2,p);
            setIntercept(fit,1);
            memcpy(XX,z_lag_1,sizeof(double)*N2);
            regress(fit,XX,z_diff,res,varcovar,alpha);
            tau = (fit->beta + 1)->value / (fit->beta + 1)->stdErr;

            scale = fit->RSS/(double)fit->df_RSS;

            p1 = 0;
            phi1_fit = reg_init(N2,p1);
            regress(phi1_fit,NULL,z_diff,res,varcovar,alpha);

            sos = phi1_fit->RSS - fit->RSS;
            dfs = (double) (phi1_fit->df_RSS - fit->df_RSS);

            phi1 = sos/dfs/scale;

            teststat[0] = tau;
            teststat[1] = phi1;
            *ltstat = 2;

            free(phi1_fit);

        } else if (!strcmp(type,"trend")) {
            p = 3;
            varcovar = (double*)malloc(sizeof(double)*p*p);
            fit = reg_init(N2,p);
            setIntercept(fit,1);
            memcpy(XX,z_lag_1,sizeof(double)*N2);
            memcpy(XX+N2,tt,sizeof(double)*N2);
            regress(fit,XX,z_diff,res,varcovar,alpha);
            tau = (fit->beta + 1)->value / (fit->beta + 1)->stdErr;

            scale = fit->RSS/(double)fit->df_RSS;

            p2 = 0;
            phi2_fit = reg_init(N2,p2);
            regress(phi2_fit,NULL,z_diff,res,varcovar,alpha);

            sos = phi2_fit->RSS - fit->RSS;
            dfs = (double) (phi2_fit->df_RSS - fit->df_RSS);

            phi2 = sos/dfs/scale;

            free(phi2_fit);

            p3 = 1;
            phi3_fit = reg_init(N2,p3);
            setIntercept(phi3_fit,1);
            regress(phi3_fit,NULL,z_diff,res,varcovar,alpha);

            sos = phi3_fit->RSS - fit->RSS;
            dfs = (double) (phi3_fit->df_RSS - fit->df_RSS);

            phi3 = sos/dfs/scale;

            teststat[0] = tau;
            teststat[1] = phi2;
            teststat[2] = phi3;
            *ltstat = 3;

            free(phi3_fit);
        }

        free(XX);
        free(varcovar);
    }

    
    if (N1 < 25){
        rowselec = 1;
    } else if (N1 < 50) {
        rowselec = 2;
    } else if (N1 < 100) {
        rowselec = 3;
    } else if (N1 < 250) {
        rowselec = 4;
    } else if (N1 < 500) {
        rowselec = 5;
    } else {
        rowselec = 6;
    }
    row1 = rowselec - 1;

    if (!strcmp(type,"none")) {
        memcpy(cval,cval_tau1+row1*3,sizeof(double)*3);
        *cvrows = 1;
        *cvcols = 3;
    } else if (!strcmp(type,"drift")) {
        memcpy(cval,cval_tau2+row1*3,sizeof(double)*3);
        memcpy(cval+3,cval_phi1+row1*3,sizeof(double)*3);
        *cvrows = 2;
        *cvcols = 3;
    } else if (!strcmp(type,"trend")) {
        memcpy(cval,cval_tau3+row1*3,sizeof(double)*3);
        memcpy(cval+3,cval_phi2+row1*3,sizeof(double)*3);
        memcpy(cval+6,cval_phi3+row1*3,sizeof(double)*3);
        *cvrows = 3;
        *cvcols = 3;
    }


    cprobs[0] = 0.01;
    cprobs[1] = 0.05;
    cprobs[2] = 0.1;

    free(x);
    free(z);
    free(tt);
    free(res);
    free(fit);
}

void ur_kpss(double *y, int N,const char* type,int lshort, int *klag, double *statistic,double *pval) {
    /*
    y - Time Series data of length N
    type - "level" or "trend"
    lshort - determines the length of lag used. lshort = 1 => lag = 4.0 * pow(((double) N)/100.0,0.25)
    else lag = 12.0 * pow(((double) N)/100.0,0.25). It can be bypassed by parameter klag
    klag - Length of the lag. Set it to > 0 instead of NULL if you wnat it to bypass lshort. Set it to NULL or negative otherwise

    Outputs

    statistic - Test statistic
    pval - p-value used to determine null hypothesis

    */
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

    if (!strcmp(type,"trend")) {
            table[0] = 0.216; table[1] = 0.176; table[2] = 0.146; table[3] = 0.119;
            p = 2;
            varcovar = (double*)malloc(sizeof(double)*p*p);

            fit = reg_init(N,p);
            regress(fit,tt,y,res,varcovar,alpha);

        } else if (!strcmp(type,"level")) {
            table[0] = 0.739; table[1] = 0.574; table[2] = 0.463; table[3] = 0.347;
            p = 1;
            varcovar = (double*)malloc(sizeof(double)*p*p);

            fit = reg_init(N,p);
            regress(fit,NULL,y,res,varcovar,alpha);
    } else {
        printf("kpss only accepts two types - level and trend. \n");
        exit(-1);
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

    if (klag != NULL && *klag > 0) {
        l = *klag;
    }

    ppsum(res,N,l,&s2);

    *statistic = eta/s2;

    arrayminmax(table,4,&ylo,&yhi);

    *pval = interpolate_linear(table,tablep,4,*statistic);

    if (klag != NULL) *klag = l;

    free(table);
    free(tt);
    free(res);
    free(csum);
    free(varcovar);
    free_reg(fit);
}

void ur_pp(double *y, int N,const char* alternative,const char* type,int lshort, int *klag, double *statistic,double *pval) {
    /*
    y - Time Series data of length N
    alternative - "stationary" or "explosive"
    type - "Z(alpha)" or "Z(t_alpha)"
    lshort - determines the length of lag used. lshort = 1 => lag = 4.0 * pow(((double) N)/100.0,0.25)
    else lag = 12.0 * pow(((double) N)/100.0,0.25). It can be bypassed by parameter klag
    klag - Length of the lag. Set it to > 0 instead of NULL if you wnat it to bypass lshort. Set it to NULL or negative otherwise

    Outputs

    statistic - Test statistic
    pval - p-value used to determine null hypothesis

    */
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

    //summary(fit);
    //anova(fit);

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

    if (klag != NULL && *klag > 0) {
        l = *klag;
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

    //printf("dx %g \n",dx);    

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

    //printf("stat %g (fit->beta+2)->value %g \n",*statistic,(fit->beta+2)->value);

    //printf("\n ss1 %g ss2 %g \n",s2,nssqr_res);

    if (klag != NULL) *klag = l;

    
    for(i = 0; i < pN;++i) {
        iter = i * tN;
        tablep[i] = interpolate_linear(nT,table+iter,tN,(double)N1);
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

    

    free_reg(fit);
    free(z);
    free(tt);
    free(XX);
    free(varcovar);
    free(res);
    free(tablep);
}

void ur_pp2(double *x, int N,const char* type,const char* model,int lshort, int *klag,double *cval,double *cprobs, double *auxstat,int *laux,double *teststat) {
    /* 
    ur_pp2
    C version of the R code by
    Bernhard Pfaff [aut, cre], Eric Zivot [ctb], Matthieu Stigler [ctb]. urca: Unit Root and Cointegration Tests for Time Series Data
    License Information : https://cran.r-project.org/web/packages/urca/urca.pdf
    */
    double *y,*y_li,*trend,*XX,*varcovar,*res,*tt,*coprods,*weights;
    int *idx;
    int N1, lmax, i,p, N12,len,iter;
    double N1d,at,my_tstat,beta_tstat,s,sum_res,sum_ymy,sum_y2,ymean,tstat;
    double myybar,myy,mty,my,sum_y,sig,lambda,lambda_prime,M,my_stat,beta_stat;
    reg_object fit;

    N1 = N - 1;
    y = &x[1];// length N1
    y_li = &x[0];// Length N1

    if (lshort == 1) {
        lmax = (int) 4.0 * pow(((double) N1)/100.0,0.25);
    } else {
        lmax = (int) 12.0 * pow(((double) N1)/100.0,0.25);
    }

    if (klag != NULL && *klag > 0) {
        lmax = *klag;
    }

    if (!strcmp(model,"trend")) {
        cval[0] = -3.9638-8.353/N1-47.44/(N1*N1);
        cval[1] = -3.4126-4.039/N1-17.83/(N1*N1);
        cval[2] = -3.1279-2.418/N1-7.58/(N1*N1);
        cprobs[0] = 0.01;
        cprobs[1] = 0.05;
        cprobs[2] = 0.1;

        //mdisplay(cval,1,3);
        trend = (double*)malloc(sizeof(double)*N1);
        tt = (double*)malloc(sizeof(double)*N1);
        N1d = (double) N1 / 2.0;

        for(i = 1; i <= N1;++i) {
            trend[i-1] = (double)i - N1d;
            tt[i-1] = (double) i;
        } 

        p = 3;
        varcovar = (double*)malloc(sizeof(double)*p*p);
        XX = (double*)malloc(sizeof(double)*2*N1);
        res = (double*)malloc(sizeof(double)*N1);

        memcpy(XX,y_li,sizeof(double)*N1);
        memcpy(XX+N1,trend,sizeof(double)*N1);

        at = 0.95;

        fit = reg_init(N1,p);

        regress(fit,XX,y,res,varcovar,at);

        /*summary(fit);
        anova(fit);
        mdisplay(res,1,N1);
        mdisplay(XX,2,N1);
        mdisplay(y,1,N1);
        */

        my_tstat = (fit->beta+0)->value/(fit->beta+0)->stdErr;
        beta_tstat = (fit->beta+2)->value/(fit->beta+2)->stdErr;

        //printf("c1 %g c2 %g \n",my_tstat,beta_tstat);

        sum_res = 0.0;
        sum_ymy = 0.0;
        sum_y2 = 0.0;
        sum_y = 0.0;
        mty = 0.0;

        ymean = mean(y,N1);
        N12 = N1 * N1;

        for(i = 0; i < N1;++i) {
            sum_res += (res[i] *res[i]);
            sum_ymy += (y[i] - ymean)*(y[i] - ymean);
            sum_y2 += (y[i] *y[i]);
            mty += (tt[i]*y[i]);
            sum_y += y[i];
        }

        s = 1.0 /(double)N1 * sum_res;
        myybar = 1.0 / (double) N12 *sum_ymy;
        myy = 1.0 / (double) N12 *sum_y2;
        mty = pow((double)N1, -2.5) * mty;
        my = pow((double)N1, -1.5) * sum_y;

        //printf("c1 %g c2 %g s %g sum_ymy %g sum_y2 %g mty %g my %g \n",my_tstat,beta_tstat,s,myybar,myy,mty,my);

        idx = (int*)malloc(sizeof(int)*lmax);
        coprods = (double*)malloc(sizeof(double)*lmax);
        weights = (double*)malloc(sizeof(double)*lmax);

        for(i = 0; i < lmax;++i) {
            idx[i] = i+1;
        }

        for(i = 0; i < lmax;++i) {
            iter = idx[i];
            len = N1 - iter;
            mmult(res+iter,res,coprods+i,1,len,1);
        }

        for(i = 0; i < lmax;++i) {
            weights[i] = 1.0 - (double) idx[i] / (double) (lmax+1);
        }

        //mdisplay(weights,1,lmax);

        mmult(weights,coprods,&sig,1,lmax,1);

        sig = s + 2.0 * sig / (double) N1;
        lambda = 0.5 * (sig - s);
        lambda_prime = lambda/sig;
        M = (1-pow((double)N1, -2.0))*myy - 12*mty*mty + 12.0*(1.0 + 1.0/(double)N1)*mty*my - (4 + 6/(double)N1 + 2/(double)(N1*N1))*my*my;
        my_stat = sqrt(s/sig)*my_tstat - lambda_prime*sqrt(sig)*my/(sqrt(M)*sqrt((M+my*my)));
        beta_stat = sqrt(s/sig)*beta_tstat - lambda_prime*sqrt(sig)*(0.5*my - mty)/(sqrt(M/12)*sqrt(myybar));
        auxstat[0] = my_stat;
        auxstat[1] = beta_stat;

        *laux = 2;

        if (!strcmp(type,"Z-tau")) {
            tstat = ((fit->beta+1)->value - 1.0)/(fit->beta+1)->stdErr;
            *teststat = sqrt(s/sig)*tstat-lambda_prime*sqrt(sig)/sqrt(M);
        } else if (!strcmp(type,"Z-alpha")) {
            *teststat = (double)N1 * ((fit->beta+1)->value - 1.0) - lambda/M;
            cval[0] = cval[1] = cval[2] = NAN;
        }

        //printf("sig %g lambda %g lambda_prime %g M %g my_stat %g beta_stat %g teststat %g \n",sig,lambda,lambda_prime,M,my_stat,beta_stat,*teststat);


        free(trend);
        free(varcovar);
        free(XX);
        free(res);
        free(tt);
        free(idx);
        free(coprods);
        free(weights);
        free_reg(fit);
    } else if (!strcmp(model,"constant")) {
        cval[0] = -3.4335-5.999/(double)N1-29.25/(double)(N1*N1);
        cval[1] = -2.8621-2.738/(double)N1-8.36/(double)(N1*N1);
        cval[2] = -2.5671-1.438/(double)N1-4.48/(double)(N1*N1);
        cprobs[0] = 0.01;
        cprobs[1] = 0.05;
        cprobs[2] = 0.1;

        //mdisplay(cval,1,3);

        p = 2;
        varcovar = (double*)malloc(sizeof(double)*p*p);
        XX = (double*)malloc(sizeof(double)*N1);
        res = (double*)malloc(sizeof(double)*N1);

        memcpy(XX,y_li,sizeof(double)*N1);

        at = 0.95;

        fit = reg_init(N1,p);

        regress(fit,XX,y,res,varcovar,at);

        //summary(fit);
        //anova(fit);
        //mdisplay(res,1,N1);
        //mdisplay(XX,2,N1);
        //mdisplay(y,1,N1);

        my_tstat = (fit->beta+0)->value/(fit->beta+0)->stdErr;

        sum_res = 0.0;
        sum_ymy = 0.0;
        sum_y2 = 0.0;
        sum_y = 0.0;

        ymean = mean(y,N1);
        N12 = N1 * N1;

        for(i = 0; i < N1;++i) {
            sum_res += (res[i] *res[i]);
            sum_ymy += (y[i] - ymean)*(y[i] - ymean);
            sum_y2 += (y[i] *y[i]);
            sum_y += y[i];
        }

        s = 1.0 /(double)N1 * sum_res;
        myybar = 1.0 / (double) N12 *sum_ymy;
        myy = 1.0 / (double) N12 *sum_y2;
        my = pow((double)N1, -1.5) * sum_y;

        //printf("c1 %g s %g sum_ymy %g sum_y2 %g my %g \n",my_tstat,s,myybar,myy,my);

        idx = (int*)malloc(sizeof(int)*lmax);
        coprods = (double*)malloc(sizeof(double)*lmax);
        weights = (double*)malloc(sizeof(double)*lmax);

        for(i = 0; i < lmax;++i) {
            idx[i] = i+1;
        }

        for(i = 0; i < lmax;++i) {
            iter = idx[i];
            len = N1 - iter;
            mmult(res+iter,res,coprods+i,1,len,1);
        }

        for(i = 0; i < lmax;++i) {
            weights[i] = 1.0 - (double) idx[i] / (double) (lmax+1);
        }

        //mdisplay(coprods,1,lmax);

        //mdisplay(weights,1,lmax);

        mmult(weights,coprods,&sig,1,lmax,1);

        sig = s + 2.0 * sig / (double) N1;
        lambda = 0.5 * (sig - s);
        lambda_prime = lambda/sig;

        my_stat = sqrt(s/sig)*my_tstat - lambda_prime*sqrt(sig)*my/(sqrt(myy)*sqrt(myybar));
        auxstat[0] = my_stat;
        *laux = 1;

        if (!strcmp(type,"Z-tau")) {
            tstat = ((fit->beta+1)->value - 1.0)/(fit->beta+1)->stdErr;
            *teststat = sqrt(s/sig)*tstat-lambda_prime*sqrt(sig)/sqrt(myybar);
        } else if (!strcmp(type,"Z-alpha")) {
            *teststat = (double)N1 * ((fit->beta+1)->value - 1.0) - lambda/myybar;
            cval[0] = cval[1] = cval[2] = NAN;
        }

        //printf("sig %g lambda %g lambda_prime %g my_stat %g teststat %g \n",sig,lambda,lambda_prime,my_stat,*teststat);
        



        free(varcovar);
        free(XX);
        free(res);
        free(idx);
        free(coprods);
        free(weights);
        free_reg(fit);
    }
}

