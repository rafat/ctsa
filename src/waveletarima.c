// SPDX-License-Identifier: BSD-3-Clause
#include "waveletarima.h"

void waveletarima(double *x, int N, char *wname, int levels,int p, int q, double *forecasts, int lengthfor) {
    wave_object obj;
	wt_object wt;
    int d, Nused, M, retval, J, i;
    double *phi, *theta;
    double *inp;

    d = 0;
    Nused = N;
    M = 1;
    retval = 0;

    J = wmaxiter(N);

    obj = wave_init(wname);
    wt = wt_init(obj, "modwt", N, J);

    modwt(wt,x);

    for(i = 0; i <= J; ++i) {
        inp = &wt->output[i*N];

    }
    

    wave_free(obj);
    wt_free(wt);
    free(phi);
    free(theta);
}