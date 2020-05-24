#include "seastest.h"

static void findMeans(double *x, int N, int stride, double *means) {
    int i, j, batches;

    batches = ceil((double) N / (double) stride);

    for(i = 0; i < stride;++i) {
        means[i] = 0.0;
    }

    for(i = 0; i < batches-1;++i) {
        for(j = 0; j < stride;++j) {
            means[j] += x[i*stride+j];
        }
    }

    for(j = 0; j < N - (batches-1)*stride;++j) {
        means[j] += x[(batches-1)*stride + j];
        means[j] /= (double) batches;
    }

    for(j = N - (batches-1)*stride; j < stride;++j) {
        means[j] /= (double) (batches - 1);
    }

}

void decompose(double *x,int N, int f,double *filter, const char *type, double *trend, int *ltrend, double *seas, int *lseas, double *random, int *lrandom) {
    int isOdd,i,j,f2,clen,batches;
    double *filt,*cout,*detrend,*seasonalmeans;


    isOdd = f%2;

    if (filter == NULL) {
        filt = (double*) malloc(sizeof(double) * f);
        for(i = 0; i < f;++i) {
            filt[i] = 1.0 / (double) f;
        }
    } else {
        filt = &filter[0];
    }

    if (!(strcmp(type,"additive") || strcmp(type,"multiplicative"))) {
        printf("Type takes only either additive or multiplicative values \n");
        exit(-1);
    }

    f2 = f / 2;
    clen = N - f + 1;

    cout = (double*) malloc(sizeof(double) * clen);

    convolve("valid","direct",x,N,filt,f,cout);

    mdisplay(cout,1,clen);

    if (!isOdd) {
        clen--;
    }

    detrend = (double*) malloc(sizeof(double) * clen);
    seasonalmeans = (double*) malloc(sizeof(double) * f);

    if (!strcmp(type,"multiplicative")) {
        for(i = f2; i < clen + f2;++i) {
            detrend[i-f2] = x[i] / cout[i-f2];
        }
    } else if (!strcmp(type,"additive")) {
        for(i = f2; i < clen + f2;++i) {
            detrend[i-f2] = x[i] - cout[i-f2];
        }
    }

    mdisplay(detrend,1,clen);

    findMeans(detrend,clen,f,seasonalmeans);

    for(i = 0; i < f-f2;++i) {
        seas[i] = seasonalmeans[f2+i];
    }

    for(i = f-f2; i < f;++i) {
        seas[i] = seasonalmeans[i-f+f2];
    }

    batches = ceil((double) clen / (double) f);

    for(i = 1; i < batches-1;++i) {
        for(j = 0; j < f;++j) {
            seas[i*f+j] =seas[j];
        }
    }
    j = 0;
    for( i = f*(batches-1); i < clen - f*(batches-1);++i) {
        seas[i] = seas[j];
        j++;
    }

    if (isOdd) {
        *lseas = clen - 1;
        *lrandom = clen - 1;
    } else {
        *lseas = clen;
        *lrandom = clen;
    }

    memcpy(trend,cout,sizeof(double)*clen);
    *ltrend = clen;

    if (!strcmp(type,"multiplicative")) {
        for(i = f2; i < clen + f2;++i) {
            random[i-f2] = x[i] / trend[i-f2];
        }
        for(i = 0; i < clen;++i) {
            random[i] = random[i] /seas[i];
        }
    } else if (!strcmp(type,"additive")) {
        for(i = f2; i < clen + f2;++i) {
            random[i-f2] = x[i] - trend[i-f2];
        }
        for(i = 0; i < clen;++i) {
            random[i] = random[i] - seas[i];
        }
    }


    if (filter == NULL) {
        free(filt);
    }

    free(cout);
    free(detrend);
    free(seasonalmeans);
}

static double calcOCSBCritVal(int f) {
    double log_f,seasonal;

    log_f = log((double)f);

    seasonal = -0.2937411*exp(-0.2850853*(log_f-0.7656451)+(-0.05983644)*((log_f-0.7656451)*(log_f-0.7656451)))-1.652202;

    return seasonal;
}

double* genLags(double *y, int N, int mlags, int *rows, int *cols) {
    int i;
    double *out;

    if (mlags <= 0) {
        out = (double*)calloc(N,sizeof(double));
        *rows = 1;
        *cols = N;
        return out;
    }

    if (mlags == 1) {
        out = (double*) malloc(sizeof(double)*N);
        memcpy(out,y,sizeof(double)*N);
        *rows = 1;
        *cols = N;
        return out;
    }

    out = (double*) malloc(sizeof(double) * (N-mlags+1)*mlags);

    for(i = 0; i < mlags;++i) {
        memcpy(out+(N-mlags+1)*i,y+mlags-1-i,sizeof(double)*(N-mlags+1));
    }

    *rows = mlags;
    *cols = N - mlags + 1;

    return out;

}

reg_object fitOCSB(double *x, int N, int f, int lag, int mlags) {
    int i,j,Nyf,Ny, ylrows, ylcols,p;
    double *y_fdiff, *y, *ylag, *mf, *res, *varcovar;
    double *z4_y,*z4_lag,*z4_preds,*z4,*inp,*z5_y,*z5_lag,*z5_preds,*z5,*XX;
    int Nz4y,Nz5y,mfrows;
    double alpha;
    reg_object fit;
    reg_object fitout;

    y_fdiff = (double*)malloc(sizeof(double)*N);
    y = (double*)malloc(sizeof(double)*N);
    double var[2] = {0,0};

    Nyf = diffs(x,N,1,f,y_fdiff);
    Ny = diff(y_fdiff,Nyf,1,y);

    //mdisplay(y_fdiff,1,Nyf);
    //mdisplay(y,1,Ny);

    ylag = genLags(y,Ny,lag,&ylrows,&ylcols);

    //mdisplay(ylag,ylrows,ylcols);

    if (mlags > -1) {
        y = &y[mlags];
        Ny -= mlags;
    }

    //mdisplay(y,1,Ny);

    mf = (double*) malloc(sizeof(double)*ylrows*Ny);
    mfrows = ylrows;

    for(i = 0; i < ylrows;++i) {
        memcpy(mf+Ny*i,ylag+ylcols*i,sizeof(double)*Ny);
    }

    mdisplay(mf,ylrows,Ny);

    p = ylrows + 1;
    varcovar = (double*)malloc(sizeof(double)*p*p);
    res = (double*)malloc(sizeof(double)*Ny);
    alpha = 0.95;
    fit = reg_init(Ny,p);

    regress(fit,mf,y,res,varcovar,alpha);

    summary(fit);
    //printf("loglik %g aic %g bic %g aicc %g \n",fit->loglik,fit->aic,fit->bic,fit->aicc);

    z4_y = &y_fdiff[lag];
    Nz4y = Nyf - lag;
    //mdisplay(z4_y,1,Nz4y);

    z4_lag = genLags(y_fdiff,Nyf,lag,&ylrows,&ylcols);

    //mdisplay(z4_lag,ylrows,ylcols);

    inp = (double*) malloc(sizeof(double)*ylcols);
    z4_preds = (double*) malloc(sizeof(double)*Nz4y);
    z4 = (double*) malloc(sizeof(double)*Nz4y);

    //printf("nz4y %d ylrows %d ylcols %d",Nz4y,ylrows,ylcols);

    for(i = 0; i < Nz4y;++i) {
        for(j = 0; j < ylrows;++j) {
            inp[j] = z4_lag[j*ylcols+i];
        }

        z4_preds[i] = fitted(fit,inp,varcovar,var);
        z4[i] = z4_y[i] - z4_preds[i];
    }

    free(inp);

    //mdisplay(z4_preds,1,Nz4y);

    z5_y = (double*)malloc(sizeof(double)*N);

    Nz5y = diff(x,N,1,z5_y);

    z5_lag = genLags(z5_y,Nz5y,lag,&ylrows,&ylcols);

    z5_y = &z5_y[lag];
    Nz5y -= lag;

    inp = (double*) malloc(sizeof(double)*ylcols);
    z5_preds = (double*) malloc(sizeof(double)*Nz5y);
    z5 = (double*) malloc(sizeof(double)*Nz5y);

    for(i = 0; i < Nz5y;++i) {
        for(j = 0; j < ylrows;++j) {
            inp[j] = z5_lag[j*ylcols+i];
        }

        z5_preds[i] = fitted(fit,inp,varcovar,var);
        z5[i] = z5_y[i] - z5_preds[i];
    }

    //mdisplay(z5_preds,1,Nz5y);

    XX = (double*) malloc(sizeof(double)*(mfrows+2) * Ny);
    memcpy(XX,mf,sizeof(double)*mfrows*Ny);
    memcpy(XX+mfrows*Ny,z4,sizeof(double)*Ny);
    memcpy(XX+(mfrows+1)*Ny,z5,sizeof(double)*Ny);

    //mdisplay(XX,mfrows+2,Ny);
    //mdisplay(y,1,Ny);

    free(varcovar);
    free(res);

    p = mfrows+2;

    fitout = reg_init(Ny,p);
    setIntercept(fitout,0);

    varcovar = (double*)malloc(sizeof(double)*p*p);
    res = (double*)malloc(sizeof(double)*Ny);

    regress(fitout,XX,y,res,varcovar,alpha);


    free(y_fdiff);
    free(y);
    free(varcovar);
    free(ylag);
    free(mf);
    free(res);
    free(fit);
    free(z4_lag);
    free(inp);
    free(z4_preds);
    free(z5_y);
    free(z5_lag);
    free(z4);
    free(z5);
    free(z5_preds);
    free(XX);

    return fitout;
}

static double getCVal(reg_object fit,const char *method) {
    double crit;

    if(!strcmp(method,"aic") || !strcmp(method,"AIC")) {
        crit = fit->aic;
    } else if(!strcmp(method,"bic") || !strcmp(method,"BIC")) {
        crit = fit->bic;
    } else if(!strcmp(method,"aicc") || !strcmp(method,"AICc") || !strcmp(method,"AICc")) {
        crit = fit->aicc;
    } else {
        printf("Only three criterions are accepted - aic, bic and aicc \n");
        exit(-1);
    }

    return crit;
}

static int checkAllNans(double *icvals,int N) {
    int i;

    for(i = 0; i < N;++i) {
        if (icvals[i] == icvals[i]) {
            return 0;
        }
    }

    return 1;
}

static int getBestIndex(double *icvals,int N) {
    int bind,i;
    double best;

    best = icvals[0] == icvals[0] ? icvals[0] : DBL_MAX;
    bind = icvals[0] == icvals[0] ? 0 : -1;

    for(i = 1; i < N;++i) {
        if (icvals[i] < best && icvals[i] == icvals[i]) {
            bind = i;
            best = icvals[i];
        }
    }

    return bind;
}

void OCSBtest(double *x, int N, int f, int mlags, const char *method) {
    int i, bestindex,allnans,maxlag;
    double stat,crit;
    double *tval,*icvals;
    reg_object fit;
    reg_object *list = (reg_object*)malloc(sizeof(reg_object)*mlags);
    reg_object crit_reg = NULL;

    icvals = (double*) malloc(sizeof(double)*mlags);

    maxlag = mlags;
    if (mlags > 0 && strcmp(method,"fixed")) {

        for(i = 1; i <= mlags;++i) {
            list[i-1] = fitOCSB(x,N,f,i,mlags);
            icvals[i-1] = getCVal(list[i-1],method);
            printf("icvals %g",icvals[i-1]);
        }

        allnans = checkAllNans(icvals,mlags);

        if (allnans == 1) {
            printf("All lag values up to 'maxlag' produced singular matrices. Use different method. \n");
            exit(-1);
        }

        bestindex = getBestIndex(icvals,mlags);
        maxlag = bestindex-1;
        crit_reg = list[bestindex];

    }

    if (maxlag <= 0) {
        fit = crit_reg;
    } else {
        fit = fitOCSB(x,N,f,maxlag,maxlag);

        if (fit->rank != fit->p && fit->rank == fit->rank) {
            if (crit_reg == NULL) {
                printf("Could not find a solution. Try a different method. \n");
                exit(-1);
            } else {
                fit = crit_reg;
            }
        }
    }

    tval = (double*) malloc(sizeof(double)*fit->p);

    for(i = 0; i < fit->p; ++i) {
		tval[i] = (fit->beta+i)->value/(fit->beta+i)->stdErr;
	}

    mdisplay(tval,1,fit->p);

    for(i = 0; i < mlags;++i) {
        free_reg(list[i]);
    }

    crit = calcOCSBCritVal(f);

    printf("crit %g \n",crit);

    free(list);
    free(icvals);
    free(tval);
}