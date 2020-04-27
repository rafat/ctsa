#include "unitroot.h"

void ur_df(double *y, int N,const char* type, int lags, const char *selectlags) {
    int lag,N1,i,j;
    double *x,*z;
    if (lags < 0) {
        printf("Lags must be >= 0 \n");
        exit(-1);
    }

    lag = lags;
    lags++;
    N1 = N - 1;

    z = (double*)malloc(sizeof(double)*N1);
    x = (double*)malloc(sizeof(double)*(N1-lags+1)*lags);

    diff(y,N,1,z);

    
    for(i = 0; i < lags;++i) {
        for(j = 0;j < N1 - lags + 1;++j) {
            x[i*(N1-lags+1)+j] = zz[lags+j-i-1];
        }
    }

    mdisplay(x,lags,N1-lags+1);
    

    free(z);
}
