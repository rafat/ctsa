// SPDX-License-Identifier: BSD-3-Clause
#include "autoutils.h"

static int runstattests(double *x, int N, const char * test,const char *mtype,double alpha) {
    int diff,klag,lshort;
    const char *alternative;
    const char *type;
    const char *model;
    const char *selectlags;
    double stat, pval;
    double cval[9] = {0,0,0,0,0,0,0,0,0};
	double cprobs[3] = {0,0,0};
	double auxstat[2] = {0,0};
    double teststat[3] = {0,0,0};
	int laux,pN,cvrows,cvcols,ltstat;

    if (!strcmp(test,"kpss")) {
        if (!strcmp(mtype,"level")) {
            type = "level";
        } else if (!strcmp(mtype,"trend")) {
            type = "trend";
        } else {
            printf("ndiffs only accepts one of two types : level and trend. \n");
            exit(-1);
        }
        
        klag = (int) 3.0 * sqrt((double)N) / 13.0;
        lshort = 1;
        ur_kpss(x,N,type,lshort,&klag,&stat,&pval);
        //printf("KPSS stat %g pval %g \n",stat,pval);
        diff = (pval < alpha) ? 1 : 0;
    } else if (!strcmp(test,"df") || !strcmp(test,"adf")) {
        if (!strcmp(mtype,"level")) {
            type = "drift";
        } else if (!strcmp(mtype,"trend")) {
            type = "trend";
        } else {
            printf("ndiffs only accepts one of two types : level and trend. \n");
            exit(-1);
        }
        selectlags = "fixed";
        klag = 1;
        pN = 3;
        //alternative = "stationary";
        //ur_df(x,N,alternative,NULL,&stat,&pval);
        ur_df2(x,N,type,&klag,selectlags,cval,&cvrows,&cvcols, cprobs, teststat,&ltstat);
        stat = teststat[0];
        pval = interpolate_linear(cval,cprobs,pN,stat);
        //printf("ADF stat %g pval %g \n",stat,pval);
        diff = (pval > alpha) ? 1 : 0;
    } else if (!strcmp(test,"pp")) {
        if (!strcmp(mtype,"level")) {
            model = "constant";
        } else if (!strcmp(mtype,"trend")) {
            model = "trend";
        } else {
            printf("ndiffs KPSS only accepts one of two types : level and trend. \n");
            exit(-1);
        }
        type = "Z-tau";
        lshort = 1;
        pN = 3;
        //ur_pp(x,N,alternative,type,lshort,NULL, &stat,&pval);
        ur_pp2(x,N,type,model,lshort,NULL,cval,cprobs,auxstat,&laux,&stat);
        pval = interpolate_linear(cval,cprobs,pN,stat);
        //printf("PP stat %g pval %g \n",stat,pval);
        diff = (pval > alpha) ? 1 : 0;
    } else {
        printf("Only three tests are allowed - kpss, df and pp \n");
        exit(-1);
    }

    return diff;
}

static int runseasonalitytests(double *x, int N,int f, const char * test) {
    int diff,mlags;
    double stat,crit;
    const char *method;

    if (!strcmp(test,"ocsb")) {
        mlags = 3;
        method = "aic";
        OCSBtest(x,N,f,mlags,method,&stat,&crit);
        diff = (stat > crit) ? 1 : 0;
    } else if (!strcmp(test,"seas")) {
        crit = 0.64;
        SHtest(x,N,&f,1,&stat);
        diff = (stat > crit) ? 1 : 0;
    }

    return diff;
}

int ndiffs(double *x, int N,double *alpha, const char *test,const char *type, int *max_d) {
    int d,max_d_,cc,NX,dodiff;
    double alpha_;
    double *y,*z;

    d = 0;

    alpha_ = alpha == NULL ? 0.05 : *alpha;
    max_d_ = max_d == NULL ? 2 : *max_d;

    // printf("%g \n",alpha_);

    if (max_d_ < 0) {
        printf("Error. Maximum Difference cannot be less than 0 \n");
        exit(-1);
    } 

    if (alpha_ < 0.01) {
        alpha_ = 0.01;
        printf("Warning : Alpha only takes values from 0.01 o 0.1. Setting alpha to 0.01. \n");
    }

    if (alpha_ > 0.1) {
        alpha_ = 0.1;
        printf("Warning : Alpha only takes values from 0.01 o 0.1. Setting alpha to 0.1. \n");
    }

    cc = checkConstant(x,N);

    if (cc == 1) return d;

    NX = N;

    dodiff = runstattests(x,N,test,type,alpha_);


    if (dodiff != dodiff) return d;

    y = (double*)malloc(sizeof(double)*NX);
    z = (double*)malloc(sizeof(double)*NX);

    memcpy(y,x,sizeof(double)*NX);

    while (dodiff && d < max_d_) {

        d++;
        //printf("d %d \n",d);

        NX = diff(y,NX,1,z);

        cc = checkConstant(z,NX);

        //printf("cc %d \n",cc);

        if (cc == 1) break;

        dodiff = runstattests(z,NX,test,type,alpha_);

        //printf("dodiff %d \n",dodiff);

        if (dodiff != dodiff) {
            d--;
            break;
        }

        memcpy(y,z,sizeof(double)*NX);
    }

    free(y);
    free(z);

    return d;
}

int nsdiffs(double *x, int N,int f,double *alpha, const char *test, int *max_D)  {
    int D,max_D_,cc,NX,dodiff;
    double alpha_;
    double *y,*z;

    D = 0;

    alpha_ = alpha == NULL ? 0.05 : *alpha;
    max_D_ = max_D == NULL ? 1 : *max_D;

    if (max_D_ < 0) {
        printf("Error. Maximum Difference cannot be less than 0 \n");
        exit(-1);
    } 

    if (alpha_ < 0.01) {
        alpha_ = 0.01;
        printf("Warning : Alpha only takes values from 0.01 o 0.1. Setting alpha to 0.01. \n");
    }

    if (alpha_ > 0.1) {
        alpha_ = 0.1;
        printf("Warning : Alpha only takes values from 0.01 o 0.1. Setting alpha to 0.1. \n");
    }

    if (!strcmp(test,"ocsb")) {
        if (alpha_ != 0.05) {
            printf("Significance levels other than 5%% are not currently supported by test='ocsb', defaulting to alpha = 0.05.");
            alpha_ = 0.05;
        }
    }

    cc = checkConstant(x,N);
    //printf("cc %d \n",cc);

    if (cc == 1) return D;

    NX = N;

    if (f <= 1) {
        printf("Warning : Only data with f > 1 can be differenced. \n");
        return 0;
    }

    dodiff = runseasonalitytests(x,N,f,test);

    y = (double*)malloc(sizeof(double)*NX);
    z = (double*)malloc(sizeof(double)*NX);

    memcpy(y,x,sizeof(double)*NX);

    while (dodiff == 1 && D < max_D_) {
        D++;
        NX = diffs(y,NX,1,f,z);

        cc = checkConstant(z,NX);

        //printf("cc %d \n",cc);

        if (cc == 1) break;

        dodiff = runseasonalitytests(z,NX,f,test);
        memcpy(y,z,sizeof(double)*NX);
    }

    free(y);
    free(z);

    return D;
}

int is_constant(double *x, int N) {
    int i, cst;
    double eps,tmp;

    eps = macheps();
    cst = 1;

    for(i = 1; i < N;++i) {
        tmp = x[i] - x[i-1];
        if (fabs(tmp) > eps) {
            return 0;
        }
    }

    return cst;
}
