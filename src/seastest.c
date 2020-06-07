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

    //mdisplay(cout,1,clen);

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

    //mdisplay(detrend,1,clen);

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

static int nextOdd(int x) {
    int y;
    if (x%2 == 0) {
        y = x+1;
    } else {
        y = x;
    }

    return y;
} 

static int degCheck(int deg) {
    if ( deg < 0 || deg > 1) {
        printf("Degree must be either 0 or 1 \n");
        exit(-1);
    }
    return deg;
}

static void cycle(int N, int f, int *cyc) {
    int i;

    for(i = 0; i < N;++i) {
        cyc[i] = i%f;
    }
}

static void applySeasonalMean(double *x,int N,int *cycle,int f) {
    int i;
    double *seasonal_means;
    int *cval;

    seasonal_means = (double*)calloc(f,sizeof(double)*f);
    cval = (int*)calloc(f,sizeof(int)*f);

    for(i = 0; i < N;++i) {
        cval[cycle[i]] += 1;
        seasonal_means[cycle[i]] += x[i];
    }

    for(i = 0; i < f;++i) {
        seasonal_means[i] /= (double) cval[i];
    }

    for(i = 0; i < N;++i) {
        x[i] = seasonal_means[cycle[i]];
    }

    free(seasonal_means);
    free(cval);
}

void stl(double *x,int N,int f, const char *s_window_type,int *s_window, int *s_degree, int *t_window, int *t_degree,int *l_window,int *l_degree,
    int *s_jump, int *t_jump, int *l_jump, int *robust,int *inner, int *outer,double *seasonal,double *trend, double *remainder) {
    
    double *rw,*work;
    int *seas_cycle;
    int i;
    int s_window_,s_degree_,t_window_,t_degree_,l_window_,l_degree_,s_jump_,t_jump_,l_jump_,inner_,outer_,periodic_,robust_;

    if (f < 2 || N < 2*f) {
        printf("Series is not periodic or has less than two periods. \n");
        exit(-1);
    }

    //Default Values

    periodic_ = 0;

    if (!strcmp(s_window_type,"period")) {
        periodic_ = 1;
        s_window_ = 10 * N + 1;
        s_degree_ = 0;
    } else {
        if (s_window == NULL) {
            printf("Error. Either set s_window_type to period or assign an integer value to s_window \n");
            exit(-1);
        } else {
            s_window_ = *s_window;
        }
    }

    if (t_window == NULL) {
        t_window_ = nextOdd(ceil( 1.5 * f / (1.0 - 1.5 / (double) s_window_)));
    }

    s_degree_ = (s_degree == NULL) ? 0 : *s_degree;

    //if (t_window == NULL) *t_window = 0;
    t_degree_ = (t_degree == NULL) ? 1 : *t_degree;
    l_window_ = (l_window == NULL) ? nextOdd(f) : *l_window;
    l_degree_ = (l_degree == NULL) ? t_degree_ : *l_degree;
    s_jump_ = (s_jump == NULL) ? ceil((double)s_window_ / 10.0) : *s_jump;
    t_jump_ = (t_jump == NULL) ? ceil((double)t_window_ / 10.0) : *t_jump;
    l_jump_ = (l_jump == NULL) ? ceil((double)l_window_ / 10.0) : *l_jump;

    robust_ = (robust == NULL) ? 0 : *robust;
    inner_ = (inner == NULL) ? (robust_ ? 1 : 2) : *inner;
    outer_ = (inner == NULL) ? (robust_ ? 15 : 0) : *outer;

    
    s_degree_ = degCheck(s_degree_);
    t_degree_ = degCheck(t_degree_);
    l_degree_ = degCheck(l_degree_);

    

    // Initialize work vectors

    rw = (double*)calloc(N,sizeof(double));
    work = (double*)calloc((N+2*f)*5,sizeof(double));
    seas_cycle = (int*)calloc(N,sizeof(int));


    stl_(x,&N,&f,&s_window_,&t_window_,&l_window_,&s_degree_,&t_degree_,&l_degree_,&s_jump_,&t_jump_,&l_jump_,&inner_,&outer_,rw,seasonal,trend,work);
	

    if (periodic_) {
        cycle(N,f,seas_cycle);
        applySeasonalMean(seasonal,N,seas_cycle,f);
    }

    for (i = 0; i < N;++i) {
        remainder[i] = x[i] - seasonal[i] - trend[i];
    }

    free(rw);
    free(work);
    free(seas_cycle);
}

void modstl(double *x, int N, int f, int *s_window,double *lambda, double *seasonal, double *trend,double *remainder) {
    int i, s_window_;
    double *y,*t;


    s_window_ = (s_window == NULL) ? 13 : *s_window;

    y = (double*)malloc(sizeof(double)*N);

    if (lambda != NULL) {
        boxcox(x,N,lambda,y);
    } else {
        memcpy(y,x,sizeof(double)*N);
    }

    if (f > 1) {
        stl(y,N,f,"null",&s_window_,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,seasonal,trend,remainder);
    } else {
        t = (double*) malloc(sizeof(double)*N);
        for(i = 0; i < N;++i) {
            t[i] = i;
        }
        supsmu(t,N,y,NULL,1,0,-1.0,trend);
        for(i = 0; i < N;++i) {
            seasonal[i] = 0.0;
            remainder[i] = y[i] - trend[i];
        }
        free(t);
    }

    free(y);
}

void mstl(double *x, int N, int *f, int *Nseas, int *s_window,double *lambda,int *iterate, double **seasonal, double *trend,double *remainder) {
    int i,j,k,s_window_,Niter,iterate_,iter,freq;
    double *y,*t,*deseas;
    double *msts,*seas;
    int *pos;

    Niter = 0;

    s_window_ = (s_window == NULL) ? 13 : *s_window;
    iterate_ = (iterate == NULL) ? 2 : *iterate;

    y = (double*)malloc(sizeof(double)*N);

    if (lambda != NULL) {
        boxcox(x,N,lambda,y);
    } else {
        memcpy(y,x,sizeof(double)*N);
    }

    msts = (double*)malloc(sizeof(double)*(*Nseas));

    for(i = 0; i < *Nseas;++i) {
        if (f[i] < N / 2) {
            msts[Niter] = (double) f[i];
            Niter++;
        }
    }

    pos = (int*)malloc(sizeof(int)*Niter);

    sort1d_ascending(msts,Niter,pos);

    if (Niter > 0) {

        seas = (double*)calloc(Niter*N,sizeof(double));
        deseas = (double*)calloc(N,sizeof(double));

        memcpy(deseas,y,sizeof(double)*N);

        for(j = 0; j < iterate_;++j) {
            for(i = 0; i < Niter;++i) {
                iter = i * N;
                for(k = 0; k < N;++k) {
                    deseas[k] += seas[iter+k];
                }
                freq = (int) msts[i];
                stl(deseas,N,freq,"null",&s_window_,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,seas+iter,trend,remainder);
                for(k = 0; k < N;++k) {
                    deseas[k] -= seas[iter+k];
                }
                f[i] = freq;
            }
        }

        for(i = 0; i < N;++i) {
            remainder[i] = deseas[i] - trend[i];
        }
        for(i = 0; i < Niter;++i) {
            iter = i * N;
            for(j = 0; j < N;++j) {
                seasonal[i][j] = seas[iter+j];
            }
        }

        *Nseas = Niter;

        free(seas);
        free(deseas);
    } else {
        t = (double*) malloc(sizeof(double)*N);
        for(i = 0; i < N;++i) {
            t[i] = i;
        }
        supsmu(t,N,y,NULL,1,0,-1.0,trend);
        for(i = 0; i < N;++i) {
            remainder[i] = y[i] - trend[i];
        }
        *Nseas = 0;
        free(t);
    }


    free(y);
    free(pos);
    free(msts);
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
    double *y2,*z5_y2;
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
        y2 = &y[mlags];
        Ny -= mlags;
    }

    //mdisplay(y,1,Ny);

    mf = (double*) malloc(sizeof(double)*ylrows*Ny);
    mfrows = ylrows;

    for(i = 0; i < ylrows;++i) {
        memcpy(mf+Ny*i,ylag+ylcols*i,sizeof(double)*Ny);
    }

    //mdisplay(mf,ylrows,Ny);

    p = ylrows + 1;
    varcovar = (double*)malloc(sizeof(double)*p*p);
    res = (double*)malloc(sizeof(double)*Ny);
    alpha = 0.95;
    fit = reg_init(Ny,p);

    regress(fit,mf,y2,res,varcovar,alpha);

    //summary(fit);
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

    z5_y2 = &z5_y[lag];
    Nz5y -= lag;

    inp = (double*) malloc(sizeof(double)*ylcols);
    z5_preds = (double*) malloc(sizeof(double)*Nz5y);
    z5 = (double*) malloc(sizeof(double)*Nz5y);

    for(i = 0; i < Nz5y;++i) {
        for(j = 0; j < ylrows;++j) {
            inp[j] = z5_lag[j*ylcols+i];
        }

        z5_preds[i] = fitted(fit,inp,varcovar,var);
        z5[i] = z5_y2[i] - z5_preds[i];
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

    regress(fitout,XX,y2,res,varcovar,alpha);

	
    free(y_fdiff);
    free(varcovar);
    free(ylag);
    free(mf);
    free(res);
    free_reg(fit);
    free(z4_lag);
    free(inp);
    free(z4_preds);
    free(y);
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

void OCSBtest(double *x, int N, int f, int mlags, const char *method,double *statistics,double *critical) {
    int i, bestindex,allnans,maxlag;
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
            //printf("icvals %g",icvals[i-1]);
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

    /* if (maxlag <= 0) {
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
    } */

    fit = fitOCSB(x,N,f,maxlag,maxlag);

    //printf("\n\n%d\n\n",fit->rank);

    if (fit->rank != fit->p && fit->rank == fit->rank) {
        if (crit_reg == NULL) {
            printf("Could not find a solution. Try a different method. \n");
            exit(-1);
        } else {
            fit = crit_reg;
        }
    }

    tval = (double*) malloc(sizeof(double)*fit->p);

    for(i = 0; i < fit->p; ++i) {
		tval[i] = (fit->beta+i)->value/(fit->beta+i)->stdErr;
	}

    //mdisplay(tval,1,fit->p);

    for(i = 0; i < mlags;++i) {
		free_reg(list[i]);
    }

    *critical = calcOCSBCritVal(f);
    *statistics = tval[fit->p-1];

    free(list);
    free(icvals);
    free(tval);
}

void SHtest(double *x, int N, int *f, int Nseas, double *season) {
    /*
    Seasonal Heuristic Test
    */
    int i,j;
    double **seasonal, *trend, *remainder, *vari;
    double vare,tmp1;

    seasonal = (double**)calloc(Nseas,sizeof(double*));
    for(i = 0; i < Nseas;++i) {
        seasonal[i] = (double*)calloc(N,sizeof(double));
    }
    trend = (double*)calloc(N,sizeof(double));
    remainder = (double*)calloc(N,sizeof(double));
    vari = (double*)calloc(N,sizeof(double));

    mstl(x, N,f,&Nseas,NULL,NULL,NULL,seasonal,trend,remainder);

    vare = var(remainder,N);

    for(j = 0; j < Nseas; ++j) {

        for(i = 0; i < N;++i) {
            vari[i] = remainder[i] + seasonal[j][i];
        }

        tmp1 = 1.0 - vare/var(vari,N);

        tmp1 = tmp1 < 1 ? tmp1 : 1;

        season[j] = tmp1 > 0 ? tmp1 : 0;

    }

    free(seasonal);
    free(trend);
    free(remainder);
    free(vari);
}

double* seasdummy(double *x, int N,int f,int *rows, int *cols) {
    int i,j,iter;
    double *tt;
    double *fmatrix,*oup;
    if ( f <= 1) {
        printf("Non-Seasonal Data \n");
        exit(-1);
    }

    tt = (double*) malloc(sizeof(double)*N);
    fmatrix = (double*)malloc(sizeof(double)*N*2*f);
    oup = (double*)malloc(sizeof(double)*N*(f-1));

    for(i = 0; i < N;++i) {
        tt[i] = (double) i + 1;
    }

    for(i = 0; i < N; ++i) {
        iter = i * 2 * f;
        for(j = 1; j <= f; ++j) {
            fmatrix[iter + 2 * j - 1 ] = sin(2 * (double) PIVAL * j * tt[i] / (double) f);
            fmatrix[iter + 2 * (j - 1) ] = cos(2 * (double) PIVAL * j * tt[i] /(double) f);
            //printf("%d %d \n",iter + 2 * j,iter + 2 * j+1);
        }
    }

    printf("PIVAL %g",PIVAL);

    //mdisplay(fmatrix,N,2*f);

    for(i = 0; i < N;++i) {
        memcpy(oup + i * (f-1),fmatrix + i * 2 * f,sizeof(double)* (f-1));
    }

    itranspose(oup,N,f-1);

    *rows = f-1;
    *cols = N;

    free(tt);
    free(fmatrix);
    return oup;
}
/* TO-DO Add Canova Hansen Test
double SDtest(double *x, int N,int f) {
    int i,j,lf,ltrunc,rows,cols,p;
    double stl,alpha;
    int *frec;
    double *R1,*res,*varcovar,*Fhat,*Fhataux;
    reg_object fit;

    if ( f <= 1) {
        printf("Non-Seasonal Data \n");
        exit(-1);
    }

    if (N <= f) {
        printf("Insufficient Data \n");
        exit(-1);
    }

    lf = (f+1) / 2;

    frec = (int*) malloc(sizeof(int)*lf);

    for(i = 0; i < lf;++i) {
        frec[i] = 1;
    }

    ltrunc = floor(f*pow((double)N/100.0,0.25));

    //printf("%d \n",ltrunc);

    R1 = seasdummy(x,N,f,&rows,&cols);

    //mdisplay(R1,rows,cols);

    // Regression

    p = rows + 1;

    varcovar = (double*)malloc(sizeof(double)*p*p);
    res = (double*)malloc(sizeof(double)*N);
    alpha = 0.95;

    fit = reg_init(N,p);

    regress(fit,R1,x,res,varcovar,alpha);

    //summary(fit);
    //anova(fit);

    //mdisplay(res,1,N);

    free(frec);
    free(R1);
    free_reg(fit);
    free(res);
    free(varcovar);
    return stl;
} 
*/
