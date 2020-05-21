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

double* genLags(double *y, int N, int maxLags, int *rows, int *cols) {
    int i;
    double *out;

    if (maxLags <= 0) {
        out = (double*)calloc(N,sizeof(double));
        *rows = 1;
        *cols = N;
        return out;
    }

    if (maxLags == 1) {
        out = (double*) malloc(sizeof(double)*N);
        memcpy(out,y,sizeof(double)*N);
        *rows = 1;
        *cols = N;
        return out;
    }

    out = (double*) malloc(sizeof(double) * (N-maxLags+1)*maxLags);

    for(i = 0; i < maxLags;++i) {
        memcpy(out+(N-maxLags+1)*i,y+maxLags-1-i,sizeof(double)*(N-maxLags+1));
    }

    *rows = maxLags;
    *cols = N - maxLags + 1;

    return out;

}