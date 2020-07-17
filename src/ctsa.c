#include "ctsa.h"

arima_object arima_init(int p, int d, int q, int N) {
	arima_object obj = NULL;
	int i,M;

	if (d > 0) {
		M = 0;
	}
	else {
		M = 1;
	}
	if (p < 0 || d < 0 || q < 0 || N <= 0) {
		printf("\n Input Values cannot be Negative. Program Exiting. \n");
		exit(-1);
	}
	//retval = 0 Input Error
	// retval = 1 Probable Success
	// retval = 4 Optimization Routine didn't converge
	// Retval = 15 Optimization Routine Encountered Inf/Nan Values
	
	obj = (arima_object)malloc(sizeof(struct arima_set) + sizeof(double) * (p+q+N-d) + sizeof(double) * (p+q+M)*(p+q+M));

	obj->p = p;
	obj->d = d;
	obj->q = q;
	obj->N = N;
	obj->Nused = N - d;
	obj->M = M;
	obj->retval = 0;
	obj->cssml = 1;

	for (i = 0; i < p + q; ++i) {
		obj->params[i] = 0.0;
	}
	obj->phi = &obj->params[0];
	obj->theta = &obj->params[p];
	obj->res = &obj->params[p + q];
	obj->vcov = &obj->params[p + q + N - d];

	obj->method = 0;// 0 - MLE, 1 - CSS, 2 - Box-Jenkins
	obj->optmethod = 5; // Default Method is 5
	obj->mean = 0.0;
	obj->var = 1.0;
	obj->lvcov = (p + q + M)*(p + q + M);
	obj->ncoeff = p + q + M;

	return obj;
}

sarima_object sarima_init(int p, int d, int q,int s, int P,int D,int Q, int N) {
	sarima_object obj = NULL;
	int i, M;

	if (d == 0 && D == 0) {
		M = 1;
	}
	else {
		M = 0;
	}

	if (p < 0 || d < 0 || q < 0 || N <= 0 || P < 0 || D < 0 || Q < 0 || s < 0) {
		printf("\n Input Values cannot be Negative. Program Exiting. \n");
		exit(-1);
	}

	//retval = 0 Input Error
	// retval = 1 Probable Success
	// retval = 4 Optimization Routine didn't converge
	// Retval = 15 Optimization Routine Encountered Inf/Nan Values

	obj = (sarima_object)malloc(sizeof(struct sarima_set) + sizeof(double)* (p + q + P + Q + N - d - s*D) + sizeof(double)* (p + q + P + Q + M)*(p + q + P + Q + M));

	obj->p = p;
	obj->d = d;
	obj->q = q;
	obj->N = N;
	obj->s = s;
	obj->P = P;
	obj->D = D;
	obj->Q = Q;
	obj->Nused = N - d - s*D;
	obj->M = M;
	obj->retval = 0;
	obj->cssml = 1;

	for (i = 0; i < p + q+P+Q; ++i) {
		obj->params[i] = 0.0;
	}
	obj->phi = &obj->params[0];
	obj->theta = &obj->params[p];
	obj->PHI = &obj->params[p + q];
	obj->THETA = &obj->params[p + q + P];
	obj->res = &obj->params[p + q+P+Q];
	obj->vcov = &obj->params[p + q + P + Q+ N - d - s*D];

	obj->method = 0;// 0 - MLE, 1 - CSS, 2 - Box-Jenkins
	obj->optmethod = 5; // Default Method is 5
	obj->mean = 0.0;
	obj->var = 1.0;
	obj->lvcov = (p + q + P + Q + M)*(p + q + P + Q + M);
	obj->ncoeff = p + q + P + Q + M;

	return obj;
}

ar_object ar_init(int method, int N) {
	ar_object obj = NULL;
	int i, ordermax,logn,p;


	if (N <= 0) {
		printf("\n Input Values cannot be Negative. Program Exiting. \n");
		exit(-1);
	}
	//retval = 0 Input Error
	// retval = 1 Probable Success
	// retval = 4 Optimization Routine didn't converge
	// Retval = 15 Optimization Routine Encountered Inf/Nan Values

	logn = (int) (10.0 * log((double)N));
	ordermax = imin(N - 1, logn);

	if (method == 2) {
		ordermax = imin(ordermax, 12); // MLE
	}

	p = ordermax;

	obj = (ar_object)malloc(sizeof(struct ar_set) + sizeof(double)* (p + N));

	obj->p = obj->order = obj->ordermax = p;
	obj->N = N;
	obj->retval = 0;
	obj->method = method;

	for (i = 0; i < p; ++i) {
		obj->params[i] = 0.0;
	}
	obj->phi = &obj->params[0];
	obj->res = &obj->params[p];

	obj->optmethod = 5;// Default Method is 5
	obj->mean = 0.0;
	obj->var = 1.0;

	return obj;
}

sarimax_object sarimax_init(int p, int d, int q,int P, int D, int Q,int s, int r,int imean, int N) {
	sarimax_object obj = NULL;
	int i,ncxreg;
	if (P + D + Q == 0) {
		s = 0;
	}

	if (imean < 0 || imean > 1) {
		printf("imean only accepts one of two values - 0 and 1 \n");
		exit(-1);
	}

	if (d + D > 0) {
		ncxreg = r;
	}
	else {
		if (imean == 0) {
			ncxreg = r;
		} else {
			ncxreg = 1 + r;
		}
	}
	if (p < 0 || d < 0 || q < 0 || N <= 0 || P < 0 || D < 0 || Q < 0) {
		printf("\n Input Values cannot be Negative. Program Exiting. \n");
		exit(-1);
	}
	/*
	Error Codes
	retval = 0 Input Error
	retval = 1 Probable Success
	retval = 4 Optimization Routine didn't converge
	retval = 7 Exogenous Variables are collinear
	retval = 10 Nonstationary AR part
	retavl = 12 Nonstationary Seasonal AR part
	retval = 15 Optimization Routine Encountered Inf/Nan Values

	*/
	
	obj = (sarimax_object)malloc(sizeof(struct sarimax_set) +
	 sizeof(double)* (p + q + P + Q + ncxreg + N - d - s*D) + sizeof(double)* (p + q + P + Q + ncxreg )*(p + q + P + Q + ncxreg ));

	obj->p = p;
	obj->d = d;
	obj->q = q;
	obj->N = N;
	obj->s = s;
	obj->P = P;
	obj->D = D;
	obj->Q = Q;
	obj->Nused = N - d - s*D;
	obj->M = ncxreg;
	obj->r = r;
	obj->retval = 0;
	obj->start = 0;
	obj->imean = imean;

	for (i = 0; i < (p + q + P + Q + ncxreg + N - d - s*D) + (p + q + P + Q + ncxreg )*(p + q + P + Q + ncxreg ); ++i) {
		obj->params[i] = 0.0;
	}
	obj->phi = &obj->params[0];
	obj->theta = &obj->params[p];
	obj->PHI = &obj->params[p + q];
	obj->THETA = &obj->params[p + q + P];
	obj->exog = &obj->params[p + q + P + Q];
	obj->res = &obj->params[p + q + P + Q + ncxreg];
	obj->vcov = &obj->params[p + q + P + Q + ncxreg + N - d - s*D];


	obj->method = 0;// 0 - MLE, 1 - CSS, 2 - Box-Jenkins
	obj->optmethod = 5; // Default Method is 5
	obj->mean = 0.0;
	obj->var = 1.0;
	obj->lvcov = (p + q + P + Q + ncxreg )*(p + q + P + Q + ncxreg );
	obj->ncoeff = p + q + P + Q + ncxreg ;

	return obj;
}

auto_arima_object auto_arima_init(int *pdqmax,int *PDQmax,int s, int r, int N) {
	auto_arima_object obj = NULL;
	int pmax,dmax,qmax,Pmax,Dmax,Qmax;
	int i,ncxreg;

	if (pdqmax == NULL) {
		pmax = 5;
		dmax = 2;
		qmax = 5;
	} else {
		pmax = pdqmax[0];
		dmax = pdqmax[1];
		qmax = pdqmax[2];
	}

	if (PDQmax == NULL) {
		Pmax = 2;
		Dmax = 1;
		Qmax = 2;
	} else {
		Pmax = PDQmax[0];
		Dmax = PDQmax[1];
		Qmax = PDQmax[2];
	}

	ncxreg = r + 2;
	if (pmax < 0 || dmax < 0 || qmax < 0 || N <= 0 || Pmax < 0 || Dmax < 0 || Qmax < 0) {
		printf("\n Input Values cannot be Negative. Program Exiting. \n");
		exit(-1);
	}

	obj = (auto_arima_object)malloc(sizeof(struct auto_arima_set) +
	 sizeof(double)* (pmax + qmax + Pmax + Qmax + ncxreg + N) + sizeof(double)* (pmax + qmax + Pmax + Qmax + ncxreg )*(pmax + qmax + Pmax + Qmax + ncxreg ));

	obj->pmax = pmax;
	obj->dmax = dmax;
	obj->qmax = qmax;
	obj->N = N;
	obj->s = s;
	obj->Pmax = Pmax;
	obj->Dmax = Dmax;
	obj->Qmax = Qmax;
	obj->M = ncxreg;
	obj->r = r;
	obj->retval = 0;
	obj->start = 0;

	for(i = 0; i < (pmax + qmax + Pmax + Qmax + ncxreg + N) + (pmax + qmax + Pmax + Qmax + ncxreg )*(pmax + qmax + Pmax + Qmax + ncxreg ); ++i ) {
		obj->params[i] = 0.0;
	}

	obj->phi = &obj->params[0];
	obj->theta = &obj->params[pmax];
	obj->PHI = &obj->params[pmax + qmax];
	obj->THETA = &obj->params[pmax + qmax + Pmax];
	obj->exog = &obj->params[pmax + qmax + Pmax + Qmax];
	obj->res = &obj->params[pmax + qmax + Pmax + Qmax + ncxreg];
	obj->vcov = &obj->params[pmax + qmax + Pmax + Qmax + ncxreg + N];


	obj->method = 0;// 0 - MLE, 1 - CSS, 2 - Box-Jenkins
	obj->optmethod = 5; // Default Method is 5
	obj->mean = 0.0;
	obj->var = 1.0;
	obj->lvcov = (pmax + qmax + Pmax + Qmax + ncxreg )*(pmax + qmax + Pmax + Qmax + ncxreg );
	obj->ncoeff = pmax + qmax + Pmax + Qmax + ncxreg ;

	// Set Default Values

	obj->stepwise = 1;
	obj->stationary = 0;
	obj->approximation = (N > 150 || s >  12) ? 1 : 0;
	//obj->approximation = 0;
	obj->Order_max = 5;
	obj->seasonal = 1;
	obj->num_models = 94;
	obj->alpha_test = 0.05;
	obj->alpha_seas = 0.05;
	obj->imean = 1;
	obj->idrift = 1;

	strcpy(obj->test,"kpss");
	strcpy(obj->seas,"seas");
	strcpy(obj->information_criteria,"aicc");
	strcpy(obj->type,"level");

	obj->p_start = 2;
	obj->q_start = 2;
	obj->P_start = 1;
	obj->Q_start = 1;


	return obj;
}

void arima_exec(arima_object obj, double *inp) {
	int p,q,d,N,M;
	double eps;

	p = obj->p;
	q = obj->q;
	d = obj->d;
	N = obj->N;

	if (obj->d > 0) {
		M = 0;
	}
	else {
		M = 1;
	}

	if (obj->method == 0) {
		obj->retval = as154(inp, obj->N, obj->optmethod, obj->p, obj->d, obj->q, obj->params, obj->params + p, &obj->mean, &obj->var, obj->params + p + q,&obj->loglik,
			obj->params + p + q + N - d,obj->cssml);
		obj->loglik = -0.5 * (obj->Nused * (2 * obj->loglik + 1.0 + log(2 * 3.14159)));
		obj->aic = -2.0 * obj->loglik + 2.0 * (obj->p + obj->q + obj->M) + 2.0;
	}
	else if (obj->method == 1) {
		obj->retval = css(inp, obj->N, obj->optmethod, obj->p, obj->d, obj->q, obj->params, obj->params + p, &obj->mean, &obj->var, obj->params + p + q, &obj->loglik,
			obj->params + p + q + N - d);
		obj->loglik = -0.5 * (obj->Nused * (2 * obj->loglik + 1.0 + log(2 * 3.14159)));
	}
	else if (obj->method == 2) {
		eps = macheps();
		obj->retval = nlalsm(inp, obj->N, obj->p, obj->d, obj->q, obj->params, obj->params + p, obj->M, &obj->mean, &obj->var, eps, obj->params + p + q + N - d, obj->params + p + q);
	}


}

void sarimax_exec(sarimax_object obj, double *inp,double *xreg)  {
	int p,q,d,N,M,P,D,Q,s,r,nd,ncxreg,cssml;
	double eps;

	p = obj->p;
	q = obj->q;
	d = obj->d;
	N = obj->N;
	P = obj->P;
	D = obj->D;
	Q = obj->Q;
	s = obj->s;
	r = obj->r;

	M = obj->M;

	if (obj->method == 0) {
		cssml = 1;
		obj->retval = as154x(inp, obj->N, xreg, obj->optmethod, obj->p, obj->d, obj->q, obj->s, obj->P, obj->D, obj->Q, obj->params, obj->params + p, obj->params + p + q, 
			obj->params + p + q + P,obj->params + p + q + P + Q, obj->r, &obj->mean, &obj->var,obj->params + p + q + P + Q + M,
			 &obj->loglik, obj->params + p + q + P + Q + M + N - d - s*D,cssml,obj->start,obj->imean);
		//mdisplay(obj->vcov,p+q+P+Q+M,p+q+P+Q+M);
		obj->loglik = -0.5 * (obj->Nused * (2 * obj->loglik + 1.0 + log(2 * 3.14159)));
		obj->aic = -2.0 * obj->loglik + 2.0 * (obj->p + obj->q + obj->P + obj->Q + obj->M) + 2.0;
	} else if (obj->method == 1) {
		cssml = 0;
		obj->retval = as154x(inp, obj->N, xreg, obj->optmethod, obj->p, obj->d, obj->q, obj->s, obj->P, obj->D, obj->Q, obj->params, obj->params + p, obj->params + p + q, 
			obj->params + p + q + P,obj->params + p + q + P + Q, obj->r, &obj->mean, &obj->var,obj->params + p + q + P + Q + M,
			 &obj->loglik, obj->params + p + q + P + Q + M + N - d - s*D,cssml,obj->start,obj->imean);
		//mdisplay(obj->vcov,p+q+P+Q+M,p+q+P+Q+M);
		obj->loglik = -0.5 * (obj->Nused * (2 * obj->loglik + 1.0 + log(2 * 3.14159)));
		obj->aic = -2.0 * obj->loglik + 2.0 * (obj->p + obj->q + obj->P + obj->Q + obj->M) + 2.0;
	} else if (obj->method == 2) {
		obj->retval = cssx(inp, obj->N, xreg, obj->optmethod, obj->p, obj->d, obj->q, obj->s, obj->P, obj->D, obj->Q, obj->params, obj->params + p, obj->params + p + q, 
			obj->params + p + q + P,obj->params + p + q + P + Q, obj->r, &obj->mean, &obj->var,&obj->loglik, obj->params + p + q + P + Q + M + N - d - s*D,obj->start,
			obj->imean);
		//mdisplay(obj->vcov,p+q+P+Q+M,p+q+P+Q+M);
		obj->loglik = -0.5 * (obj->Nused * (2 * obj->loglik + 1.0 + log(2 * 3.14159)));
	} else {
		printf("Only three methods are supported : 0 , 1 and 2 , where 0 is CSS-MLE ,1 is MLE and 2 is CSS \n");
		exit(-1);
	}
}

void auto_arima_exec(auto_arima_object obj, double *inp,double *xreg) {
	aa_ret_object fit;
	int p,d,q,P,D,Q,s,r,N,M;
	int i,imean,idrift, retval, iter;
	int order[3];
	int seasonal[3];
	int start[4];

	order[0] = obj->pmax;
	order[1] = obj->dmax;
	order[2] = obj->qmax;
	seasonal[0] = obj->Pmax;
	seasonal[1] = obj->Dmax;
	seasonal[2] = obj->Qmax;
	start[0] = obj->p_start;
	start[1] = obj->q_start;
	start[2] = obj->P_start;
	start[3] = obj->Q_start;

	fit = auto_arima1(inp,obj->N,order,seasonal,&obj->Order_max,obj->s,NULL,NULL,start,&obj->stationary,&obj->seasonal, obj->information_criteria,
	 &obj->stepwise,&obj->num_models,&obj->approximation,&obj->method,xreg,obj->r,obj->test,obj->type, &obj->alpha_test,obj->seas, &obj->alpha_seas,
	 &obj->idrift, &obj->imean, NULL);

	imean = obj->imean;
	idrift = obj->idrift;

	if (fit->otype == 2) {
		p = fit->Arima->sarimax->p;
		d = fit->Arima->sarimax->d;
		q = fit->Arima->sarimax->q;
		P = fit->Arima->sarimax->P;
		D = fit->Arima->sarimax->D;
		Q = fit->Arima->sarimax->Q;
		r = fit->Arima->sarimax->r;
		s = fit->Arima->sarimax->s;
		M = fit->Arima->sarimax->M;
		retval = fit->Arima->sarimax->retval;
		obj->mean = fit->Arima->sarimax->mean;
	    obj->var = fit->Arima->sarimax->var;
		obj->loglik = fit->Arima->sarimax->loglik;

		iter = (p + q + P + Q + M + N - d - s*D) + (p + q + P + Q + M )*(p + q + P + Q + M );

		for(i = 0; i < iter; ++i) {
			obj->params[i] = fit->Arima->sarimax->params[i];
		}

		obj->sigma2 = fit->Arima->sigma2;
		obj->aic = fit->Arima->aic;
		obj->bic = fit->Arima->bic;
		obj->aicc = fit->Arima->aicc;
		obj->sigma2 = fit->Arima->sigma2;
		
	} else if (fit->otype == 1) {
		p = fit->myarima->sarimax->p;
		d = fit->myarima->sarimax->d;
		q = fit->myarima->sarimax->q;
		P = fit->myarima->sarimax->P;
		D = fit->myarima->sarimax->D;
		Q = fit->myarima->sarimax->Q;
		r = fit->myarima->sarimax->r;
		s = fit->myarima->sarimax->s;
		M = fit->myarima->sarimax->M;
		retval = fit->myarima->sarimax->retval;
		obj->mean = fit->myarima->sarimax->mean;
	    obj->var = fit->myarima->sarimax->var;
		obj->loglik = fit->myarima->sarimax->loglik;

		iter = (p + q + P + Q + M + N - d - s*D) + (p + q + P + Q + M )*(p + q + P + Q + M );

		for(i = 0; i < iter; ++i) {
			obj->params[i] = fit->myarima->sarimax->params[i];
		}
		obj->sigma2 = fit->myarima->sigma2;
		obj->aic = fit->myarima->aic;
		obj->bic = fit->myarima->bic;
		obj->aicc = fit->myarima->aicc;
		obj->sigma2 = fit->myarima->sigma2;
	}


	obj->p = p;
	obj->d = d;
	obj->q = q;
	obj->N = N;
	obj->s = s;
	obj->P = P;
	obj->D = D;
	obj->Q = Q;
	obj->Nused = N - d - s*D;
	obj->M = M;
	obj->r = r;
	obj->retval = retval;
	

	obj->phi = &obj->params[0];
	obj->theta = &obj->params[p];
	obj->PHI = &obj->params[p + q];
	obj->THETA = &obj->params[p + q + P];
	obj->exog = &obj->params[p + q + P + Q];
	obj->res = &obj->params[p + q + P + Q + M];
	obj->vcov = &obj->params[p + q + P + Q + M + N - d - s*D];

	obj->lvcov = (p + q + P + Q + M )*(p + q + P + Q + M );
	obj->ncoeff = p + q + P + Q + M ;

	aa_ret_summary(fit);

	//aa_ret_free(fit);
}

void sarima_exec(sarima_object obj, double *inp) {
	int p, q, d, s,P,D,Q,N, M,cssml,trparams;
	double eps;

	p = obj->p;
	q = obj->q;
	d = obj->d;
	N = obj->N;
	P = obj->P;
	Q = obj->Q;
	D = obj->D;
	s = obj->s;

	M = obj->M;

	if (obj->method == 0) {
		obj->retval = as154_seas(inp, obj->N, obj->optmethod, obj->p, obj->d, obj->q, obj->s, obj->P, obj->D, obj->Q, obj->params, obj->params + p, obj->params + p + q, 
			obj->params + p + q + P, &obj->mean, &obj->var, &obj->loglik, obj->params + p + q + P + Q + N - d - s*D,obj->cssml);
		obj->loglik = -0.5 * (obj->Nused * (2 * obj->loglik + 1.0 + log(2 * 3.14159)));
		obj->aic = -2.0 * obj->loglik + 2.0 * (obj->p + obj->q + obj->P + obj->Q + obj->M) + 2.0;
	}
	else if (obj->method == 1) {
		obj->retval = css_seas(inp, obj->N, obj->optmethod, obj->p, obj->d, obj->q, obj->s, obj->P, obj->D, obj->Q, obj->params, obj->params + p, obj->params + p + q,
			obj->params + p + q + P, &obj->mean, &obj->var, &obj->loglik, obj->params + p + q + P + Q + N - d - s*D);
		obj->loglik = -0.5 * (obj->Nused * (2 * obj->loglik + 1.0 + log(2 * 3.14159)));
	}
	else if (obj->method == 2) {
		eps = macheps();
		obj->retval = nlalsms(inp, obj->N, obj->p, obj->d, obj->q, obj->params, obj->params + p, obj->s, obj->P, obj->D, obj->Q, obj->params + p + q, obj->params + p + q + P, obj->M, &obj->mean, &obj->var, eps,
			obj->params + p + q + P + Q + N - d - s*D, obj->params + p + q + P + Q);
	}

}

sarimax_wrapper_object sarimax_wrapper(sarimax_wrapper_object model,double *y, int N,int *order, int *seasonal, double *xreg, int r, int idrift,int mean,
	double *lambda, int biasadj,int method) {
	/*
	Init here. { sarimax_object,drift}
	*/
	sarimax_wrapper_object obj = NULL;
	int i, p,d,q,P,D,Q,s,ncoeff,drift,rr;
	double *x, *origx,*xreg2;
	double rsum;

	x = (double*) malloc(sizeof(double)*N);
	origx = (double*) malloc(sizeof(double)*N);
	obj = (sarimax_wrapper_object) malloc (sizeof(struct sarimax_wrapper_set));

	memcpy(x,y,sizeof(double)*N);
	memcpy(origx,y,sizeof(double)*N);

	if (lambda != NULL) boxcox_eval(y,N,*lambda,x);

	if (seasonal == NULL) {
		P = D = Q = s = 0;
	}


	if (model == NULL) {
		p = order[0];
		d = order[1];
		q = order[2];
		if (seasonal != NULL) {
			P = seasonal[0];
			D = seasonal[1];
			Q = seasonal[2];
			s = seasonal[3];
		}

		if (idrift == 1) {
			drift = 1;
		} else {
			drift = 0;
		}

		if (d + D > 1 && idrift == 1) {
			drift = 0;
		}

		printf("drift %d \n",drift);

		rr = r;

		
		if (drift == 1) {
			obj->idrift = 1;
			rr++;
			xreg2 = (double*)malloc(sizeof(double)*N*rr);
			for(i = 0; i < N;++i) {
				xreg2[i] = (double) (i+1);
			}
			memcpy(xreg2+N,xreg,sizeof(double)*N*(rr-1));
		} else {
			xreg2 = (double*)malloc(sizeof(double)*N*rr);
			memcpy(xreg2,xreg,sizeof(double)*N*rr);
		}
		
		obj->sarimax = sarimax_init(p,d,q,P,D,Q,s,rr,mean,N);

		sarimax_exec(obj->sarimax,x,xreg2);
		obj->aic = obj->sarimax->aic;
		ncoeff = (obj->sarimax->p + obj->sarimax->q + obj->sarimax->P + obj->sarimax->Q + obj->sarimax->M) + 1;
		obj->aicc = obj->aic + 2 * ncoeff * ((double)obj->sarimax->Nused / (double) (obj->sarimax->Nused - ncoeff - 1) - 1.0);
		obj->bic = obj->aic + ncoeff * (log((double)obj->sarimax->Nused) - 2.0);
		
		rsum = 0.0;

		for(i = 0; i < obj->sarimax->Nused;++i) {
			rsum += obj->sarimax->res[i]*obj->sarimax->res[i];
		}

		obj->sigma2 = rsum / (double) (obj->sarimax->Nused - ncoeff + 1);
	}

	free(x);
	free(origx);

	return obj;
}

void arima2(sarimax_wrapper_object model,double *x, int N,int drift,double *xreg, int r,int method) {
	int i;
	double sigma2;

	sigma2 = model->sigma2;

}

myarima_object myarima(double *x, int N, int *order, int *seasonal, int constant, const char* ic, int trace, int approx,
	double offset, double *xreg, int r, int *method) {

	myarima_object fit = NULL;
	int m,p,d,q,P,D,Q,s,diffs,i,j,imean,ncoeff,isum,imx,rr;
	int use_season,rmethod,last_nonzero,ip,iq,ir,retval;
	double rsum,minroot,tol,temp;
	double *xreg2,*phi,*theta,*zeror,*zeroi;
	int *K;

	fit = (myarima_object) malloc (sizeof(struct myarima_set));

	tol = 1e-08;

	if (order) {
		p = order[0];
		d = order[1];
		q = order[2];
	} else {
		p = 0;
		d = 0;
		q = 0;
	}

	if (seasonal) {
		P = seasonal[0];
		D = seasonal[1];
		Q = seasonal[2];
		s = seasonal[3];
		if ( s > 0) {
			m = s;
		} else {
			m = 1;
		}
	} else {
		P = 0;
		D = 0;
		Q = 0;
		s = 0;
		m = 1;
	}

	diffs = d + D;

	if ( P + D + Q > 0 && s > 0) {
		use_season = 1;
	} else {
		use_season = 0;
	}

	if (method == NULL) {
		if (approx) {
			rmethod = 2;
		} else {
			rmethod = 0;
		}
	} else {
		rmethod = *method;
	}

	rr = r;

	if (diffs == 1 && constant == 1) {
		imean = 1;
		rr++;
		xreg2 = (double*)malloc(sizeof(double)*N*rr);
		for(i = 0; i < N;++i) {
			xreg2[i] = (double) (i+1);
		}
		memcpy(xreg2+N,xreg,sizeof(double)*N*(rr-1));
		fit->sarimax = sarimax_init(p,d,q,P,D,Q,s,rr,imean,N);
		sarimax_exec(fit->sarimax,x,xreg2);
		free(xreg2);
		//sarimax_summary(fit->sarimax);
	} else {
		imean = constant;
		if (rr > 0) {
			xreg2 = (double*)malloc(sizeof(double)*N*rr);
			memcpy(xreg2,xreg,sizeof(double)*N*rr);
			fit->sarimax = sarimax_init(p,d,q,P,D,Q,s,rr,imean,N);
			sarimax_setMethod(fit->sarimax,rmethod);
			sarimax_exec(fit->sarimax,x,xreg2);
			free(xreg2);
		} else {
			printf("p: %d d: %d q: %d P: %d D: %d Q: %d \n",p,d,q,P,D,Q);
			fit->sarimax = sarimax_init(p,d,q,P,D,Q,s,0,imean,N);
			sarimax_setMethod(fit->sarimax,rmethod);
			sarimax_exec(fit->sarimax,x,NULL);
			//sarimax_summary(fit->sarimax);
		}
		
	}

	if (fit->sarimax->retval == 1) {
		ncoeff = (fit->sarimax->p + fit->sarimax->q + fit->sarimax->P + fit->sarimax->Q + fit->sarimax->M) + 1;

		if (rmethod == 2) {
			fit->sarimax->aic = offset + fit->sarimax->Nused * log(fit->sarimax->var) + 2 * ncoeff;// sigma2 
		}

		if (fit->sarimax->aic) {
			fit->aic = fit->sarimax->aic;
			fit->bic = fit->aic + ncoeff * (log((double)fit->sarimax->Nused) - 2.0);
			fit->aicc = fit->aic + 2 * ncoeff * (ncoeff + 1) / (double) (fit->sarimax->Nused - ncoeff - 1);
			if (!strcmp(ic,"aic")) {
				fit->ic = fit->aic;
			} else if (!strcmp(ic,"bic")) {
				fit->ic = fit->bic;
			} else if (!strcmp(ic,"aicc")) {
				fit->ic = fit->aicc;
			}
		} else {
			fit->ic = fit->aic = fit->bic = fit->aicc = DBL_MAX;
		}

		rsum = 0.0;

		for(i = 0; i < fit->sarimax->Nused;++i) {
			rsum += fit->sarimax->res[i]*fit->sarimax->res[i];
		}

		fit->sigma2 = rsum / (double) (fit->sarimax->Nused - ncoeff + 1);

		minroot = 2.0;

		ip = p + s*P + 1;
		iq = q + s*Q + 1;

		phi = (double*)calloc(ip,sizeof(double));
		theta = (double*)calloc(iq,sizeof(double));

		phi[0] = theta[0] = 1.0;

		for(i = 0; i < p;++i) {
			phi[i+1] = fit->sarimax->phi[i];
		}

		for(i = 0; i < q; ++i) {
			theta[i+1] = -1.0 * fit->sarimax->theta[i];
		}

		for (j = 0; j < P; ++j) {
			phi[(j + 1)*s] += fit->sarimax->PHI[j];
			for (i = 0; i < p; ++i) {
				phi[(j + 1)*s + i + 1] -= fit->sarimax->phi[i] * fit->sarimax->PHI[j];
			}
		}

		for (j = 0; j < Q; ++j) {
			theta[(j + 1)*s] -= fit->sarimax->THETA[j];
			for (i = 0; i < q; ++i) {
				theta[(j + 1)*s + i + 1] += fit->sarimax->theta[i] * fit->sarimax->THETA[j];
			}
		}

		ir = (ip > iq) ? ip : iq;

		K = (int*)malloc(sizeof(int)*ir);
		zeror = (double*)malloc(sizeof(double)*ir);
		zeroi = (double*)malloc(sizeof(double)*ir);

		if (p + P > 0) {
			isum = 0;
			imx = -1;
			for(i = 0; i < ip - 1; ++i) {
				K[i] = (fabs(phi[i+1]) > tol) ? 1 : 0; 
				if (K[i] == 1) {
					isum += 1;
					imx = i;
				}
			}

			last_nonzero = imx + 1;

			if (last_nonzero > 0) {
				for(i = 0; i < last_nonzero; ++i) {
					phi[i+1] *= -1.0;
				}

				retval = polyroot(phi,last_nonzero,zeror,zeroi);

				if (retval == 0) {
					for(i = 0; i < last_nonzero;++i) {
						temp = sqrt(zeror[i]*zeror[i] + zeroi[i]*zeroi[i]);
						if (temp < minroot) minroot = temp;
					}
				} else {
					fit->ic = DBL_MAX;
				}

			}

		}

		if (q + Q > 1 && fit->ic < DBL_MAX) {
			isum = 0;
			imx = -1;
			for(i = 0; i < iq - 1; ++i) {
				K[i] = (fabs(theta[i+1]) > tol) ? 1 : 0; 
				if (K[i] == 1) {
					isum += 1;
					imx = i;
				}
			}

			last_nonzero = imx + 1;

			if (last_nonzero > 0) {

				retval = polyroot(theta,last_nonzero,zeror,zeroi);

				if (retval == 0) {
					for(i = 0; i < last_nonzero;++i) {
						temp = sqrt(zeror[i]*zeror[i] + zeroi[i]*zeroi[i]);
						if (temp < minroot) minroot = temp;
					}
				} else {
					fit->ic = DBL_MAX;
				}

			}

			
		}

		if (minroot < 1.01) {
			fit->ic = DBL_MAX;
		}



		free(phi);
		free(theta);
		free(K);
		free(zeror);
		free(zeroi);
	} else {
		fit->ic = DBL_MAX;
	}

	return fit;
}

myarima_object search_arima(double *x, int N,int d, int D, int p_max, int q_max, int P_max, int Q_max, int Order_max, int stationary,int s, const char *ic,
	int approximation, double *xreg, int r, double offset,int allowdrift, int allowmean, int method) {

	int idrift, imean, maxK, i, j, I, J,K;
	double best_ic,bestK,iapprox;
	myarima_object bestfit;
	myarima_object fit;
	int order[3] = {0,0,0};
	int seasonal[4] = {0,0,0,0};
	int bestorder[3] = {0,0,0};
	int bestseasonal[4] = {0,0,0,0};
	int trace = 0;

	
	idrift = allowdrift && (d + D == 1);

	imean = allowmean && (d + D == 0);

	maxK = (idrift || imean);

	//serial implementation

	best_ic = DBL_MAX;

	for(i = 0; i <= p_max; ++i) {
		for(j = 0;j <= q_max;++j) {
			for(I = 0; I <= P_max;++I) {
				for(J = 0; J <= Q_max; ++J) {
					if (i+j+I+J <= Order_max) {
						for(K = 0; K <= maxK; ++K) {
							order[0] = i;
							order[1] = d;
							order[2] = j;
							
							if (s == 0) {
								seasonal[0] = seasonal[1] = seasonal[2] = seasonal[3] = 0;
							} else if (I == J == D == 0) {
								seasonal[3] = 0;
							} else {
								seasonal[0] = I;
								seasonal[1] = D;
								seasonal[2] = J;
								seasonal[3] = s;
							}

							fit = myarima(x,N,order,seasonal, K, ic, trace, approximation, offset,xreg, r, &method);
							myarima_summary(fit);

							if (best_ic > fit->ic) {
								best_ic = fit->ic;
								memcpy(bestorder,order,sizeof(int)*3);
								memcpy(bestseasonal,seasonal,sizeof(int)*4);
								bestK = K;
							}

							myarima_free(fit);

						}
					}
				}
			}
		}
	}

	printf("best_ic %g \n",best_ic);

	bestfit = myarima(x,N,bestorder,bestseasonal, bestK, ic, trace, approximation, offset,xreg, r, &method);

	if (approximation) {
		if (bestfit->ic == DBL_MAX) {
			iapprox = 0;
			myarima_free(bestfit);

			bestfit = search_arima(x,N,d,D,p_max,q_max,P_max,Q_max,Order_max,stationary,s,ic,iapprox,xreg,r,offset,allowdrift,allowmean,method);

		}
	}


	return bestfit;
}

static void set_results(double *results, int row, int p, int d, int q, int P, int D, int Q, int constant, double ic) {
	int index;

	index = 8 * (row-1);

	results[index] = (double) p;
	results[index+1] = (double) d;
	results[index+2] = (double) q;
	results[index+3] = (double) P;
	results[index+4] = (double) D;
	results[index+5] = (double) Q;
	results[index+6] = (double) constant;
	results[index+7] = ic;
}


static int newmodel(int p, int d, int q, int P, int D, int Q, int constant, double *result, int n) {
	int i;

	for(i = 0; i < n;++i) {
		if ((double)p == result[0] && (double)d == result[1] && (double)q == result[2] && (double)P == result[3] &&
		 (double)D == result[4] && (double)Q == result[5] && (double)constant == result[6]) {
			return 0;
		}
	}

	return 1;
}

aa_ret_object auto_arima1(double *y, int N, int *ordermax, int *seasonalmax,int *maxcoeff, int s,int *DD, int *dd, int *start, int *stationary, int *seasonal, 
	const char *ic, int *stepwise, int *nmodels,int *approximation,int *method,double *xreg, int r, const char *test,const char *type, double *test_alpha, 
	const char *seas, double *seas_alpha, int *allowdrift, int *allowmean, double *lambda) {
	
	
	aa_ret_object fit = NULL;
	int p_max,d_max,q_max,P_max, D_max,Q_max,p,q,P,Q,trace;
	int p_start, q_start, P_start, Q_start,d,D,Nd,Ndd,amethod, iapprox,constant,Order_max;
	int istationary, iseasonal, istepwise, models, idrift, imean, m, N3, N3m,rnk, i, is1,k;
	double *x,*xx,*varcovar,*diffxreg,*dx,*diffdxreg,*diffdx,*results,*icvector;
	reg_object reg;
	sarimax_object approxfit;
	myarima_object bestfit;
	int order[3] = {0,0,0};
	int seasonalorder[4] = {0,0,0,0};
	int bestorder[3] = {0,0,0};
	int bestseasonalorder[4] = {0,0,0,0};
	int biasadj = 0,bestconstant,startk,ntconstant;
	double offset, best_ic;
	int *icindex;

	printf("DEBUG 1. \n");

	fit = (aa_ret_object) malloc (sizeof(struct aa_ret_set));

	fit->otype = 0;

	if (ordermax == NULL) {
		p_max = 5;
		q_max = 5;
		d_max = 2;
	} else {
		p_max = ordermax[0];
		d_max = ordermax[1];
		q_max = ordermax[2];
	}

	if (seasonalmax == NULL) {
		P_max = 2;
		Q_max = 2;
		D_max = 1;
	} else {
		P_max = seasonalmax[0];
		D_max = seasonalmax[1];
		Q_max = seasonalmax[2];
	}

	if (maxcoeff == NULL) {
		Order_max = 5;
	} else {
		Order_max = *maxcoeff;
	}

	if (start == NULL) {
		p_start = 2;
		q_start = 2;
		P_start = 1;
		Q_start = 1;
	} else {
		p_start = start[0];
		q_start = start[1];
		P_start = start[2];
		Q_start = start[3];
	}

	if (stationary == NULL) {
		istationary = 0;
	} else {
		istationary = *stationary;
	}

	if (seasonal == NULL) {
		iseasonal = 1;
	}

	if (stepwise == NULL) {
		istepwise = 1;
	} else {
		istepwise = *stepwise;
	}

	if (nmodels == NULL) {
		models = 94;
	} else {
		models = *nmodels;
	}

	if (allowmean == NULL) {
		imean = 1;
	} else {
		imean = *allowmean;
	}

	if (allowdrift == NULL) {
		idrift = 0;
	} else {
		idrift = *allowdrift;
	}

	if (method == NULL) {
		amethod = 0;
	} else {
		amethod = *method;
	}


	x = (double*) malloc(sizeof(double)*N);

	memcpy(x,y,sizeof(double)*N);

	mdisplay(x,1,N);

	printf("DEBUG 2. \n");

	if (s <= 1) {
		m = 1;
	} else {
		m = s;
	}

	if (approximation == NULL) {
		if (s > 12 || N > 150) {
			iapprox = 1;
		} else {
			iapprox = 0;
		}
	} else {
		iapprox = *approximation;
	}

	N3 = N / 3;

	p_max = (p_max < N3) ? p_max : N3;
	q_max = (q_max < N3) ? q_max : N3;

	N3 = N3 / m;

	P_max = (P_max < N3m) ? P_max : N3m;
	Q_max = (Q_max < N3m) ? Q_max : N3m;

	if ( N <= 3) ic = "aic";

	if (lambda != NULL) boxcox_eval(y,N,*lambda,x);

	mdisplay(x,1,N);

	xx = (double*) malloc(sizeof(double)*N);
	dx = (double*) malloc(sizeof(double)*N);
	diffxreg = (double*) malloc(sizeof(double)*N*r);
	diffdxreg = (double*) malloc(sizeof(double)*N*r);
	diffdx = (double*) malloc(sizeof(double)*N);

	memcpy(xx,x,sizeof(double)*N);

	if (xreg) {
		// Collinearity check

		rnk = rank(xreg,N,r);

		if (rnk < r) {
			free(xx);
			free(x);
			printf("Exogenous Variables are collinear. \n");
			exit(-1);
		}

		varcovar = (double*)malloc(sizeof(double)*r*r);

		reg = reg_init(N,r);

		//setIntercept(reg,0);
		regress(reg,xreg,x,xx,varcovar,0.95);

		free(varcovar);
	}

	printf("DEBUG 3. \n");

	D = DD == NULL ? -1 : *DD;

	if (istationary) {
		d = D = 0;
	}

	if (m == 1) {
		D = P_max = Q_max = 0;
	} else if (D == -1) {
		D = nsdiffs(xx,N,s,seas_alpha,seas,&D_max);

		if (D > 0 && xreg) {
			is1 = 0;
			Nd = N - s * D;
			for(i = 0; i < r;++i) {
				diffs(xreg+i*N,N,D,s,diffxreg+i*Nd);
				is1 += is_constant(diffxreg+i*Nd,Nd);
			}
			
			
			if (is1 > 0) {
				D = D - 1;
			}
		}
	}

	printf("D %d \n",D);

	Nd = N - s * D;

	if (D > 0) {
		Nd = diffs(xx,N,D,s,dx);
	} else {
		memcpy(dx,xx,sizeof(double)*N);
		Nd = N;
	}

	if (xreg) {
		if (D > 0) {
			diffs(xreg,N,D,s,diffxreg);
		} else {
			memcpy(diffxreg,xreg,N*r);
		}
	}

	d = dd == NULL ? -1 : *dd;

	printf("d %d Nd %d \n",d,Nd);

	if (d == -1) {
		printf("alpha %g test %s type %s d_max %d \n",test_alpha,test,type,d_max);
		d = ndiffs(dx,Nd,test_alpha,test,type,&d_max);

		Ndd = Nd - d;

		if (d > 0 && xreg) {
			is1 = 0;
			for(i = 0; i < r;++i) {
				diff(diffxreg+i*Nd,Nd,d,diffdxreg+i*Ndd);
				is1 += is_constant(diffdxreg+i*Ndd,Ndd);
			}
			if (is1 > 0) {
				d = d - 1;
			}
		}
	}

	printf("d %d \n",d);

	Ndd = Nd - d;

	if (d > 0) {
		Ndd = diff(dx,Nd,d,diffdx);
	} else {
		Ndd = Nd;
		memcpy(diffdx,dx,sizeof(double)*Ndd);
	}

	printf("DEBUG 4. \n");

	if (Ndd == 0) {
		free(x);
		free(xx);
		free(dx);
		free(diffxreg);
		free(diffdxreg);
		free(diffdx);
		fit->Arima = NULL;
		fit->myarima = NULL;
		printf("Warning. Data length is 0 after differencing. Exiting. \n");
		return fit;
	} else if (is_constant(diffdx,Ndd)) {
		if (xreg == NULL) {
			if (D > 0 && d == 0) {
				seasonalorder[0] = seasonalorder[2] = 0;
				seasonalorder[1] = D;
				seasonalorder[3] = m;
				order[0] = order[1] = order[2] = 0;
				imean = 1;
				idrift = 1;
				fit->Arima = sarimax_wrapper(NULL,x,N,order,seasonal,NULL,r,idrift,imean,NULL,biasadj,amethod);
				fit->myarima = NULL;
				fit->otype = 2;
			}
			else if ( D > 0 && d > 0) {
				seasonalorder[0] = seasonalorder[2] = 0;
				seasonalorder[1] = D;
				seasonalorder[3] = m;
				order[0] = order[2] = 0;
				order[1] = d;
				idrift = 0;
				fit->Arima = sarimax_wrapper(NULL,x,N,order,seasonal,NULL,r,idrift,imean,NULL,biasadj,amethod);
				fit->myarima = NULL;
				fit->otype = 2;
			} else if (d == 2) {
				seasonalorder[0] = seasonalorder[2] = 0;
				seasonalorder[1] = 0;
				seasonalorder[3] = 0;
				order[0] = order[2] = 0;
				order[1] = d;
				idrift = 0;
				fit->Arima = sarimax_wrapper(NULL,x,N,order,seasonal,NULL,r,idrift,imean,NULL,biasadj,amethod);
				fit->myarima = NULL;
				fit->otype = 2;
			} else if (d < 2) {
				seasonalorder[0] = seasonalorder[2] = 0;
				seasonalorder[1] = 0;
				seasonalorder[3] = 0;
				order[0] = order[2] = 0;
				order[1] = d;
				idrift = 1;
				fit->Arima = sarimax_wrapper(NULL,x,N,order,seasonal,NULL,r,idrift,imean,NULL,biasadj,amethod);
				fit->myarima = NULL;
				fit->otype = 2;
			} else {
				printf("data is not sitable for ARIMA modelling. \n");
				fit->Arima = NULL;
				fit->myarima = NULL;
			}
		} else {
			if (D > 0) {
				seasonalorder[0] = seasonalorder[2] = 0;
				seasonalorder[1] = D;
				seasonalorder[3] = m;
				order[0] = order[2] = 0;
				order[1] = d;
				fit->Arima = sarimax_wrapper(NULL,x,N,order,seasonal,xreg,r,idrift,imean,NULL,biasadj,amethod);
				fit->myarima = NULL;
				fit->otype = 2;
			} else {
				seasonalorder[0] = seasonalorder[2] = 0;
				seasonalorder[1] = 0;
				seasonalorder[3] = 0;
				order[0] = order[2] = 0;
				order[1] = d;
				fit->Arima = sarimax_wrapper(NULL,x,N,order,seasonal,xreg,r,idrift,imean,NULL,biasadj,amethod);
				fit->myarima = NULL;
				fit->otype = 2;
			}

		}

		free(x);
		free(xx);
		free(dx);
		free(diffxreg);
		free(diffdxreg);
		free(diffdx);
		return fit;

	}


	if (m > 1) {
		if (P_max > 0) {
			p_max = (p_max < m-1) ? p_max : m - 1;
		}
		if (Q_max > 0) {
			q_max = (q_max < m-1) ? q_max : m - 1;
		}
	}

	printf("iapprox %d \n",iapprox);

	if (iapprox) {
		approxfit = sarimax_init(0,d,0,0,D,0,s,r,imean,N);
		mdisplay(x,1,N);
		if (r == 0) {
			sarimax_exec(approxfit,x,NULL);
		} else {
			sarimax_exec(approxfit,x,xreg);
		}

		if (approxfit->retval == 1) {
			offset = -2.0 * approxfit->loglik - (double) N * log(approxfit->var);
		} else {
			offset = 0;
		}

		sarimax_summary(approxfit);

		sarimax_free(approxfit);
	} else {
		offset = 0;
	}

	printf("offset %g \n",offset);

	idrift = idrift && (d + D == 1);
	imean = imean && (d + D == 0);

	constant =  (idrift || imean);

	if (!istepwise) {
		fit->otype = 1;
		fit->Arima = NULL;

		imean = 1;
	
		fit->myarima = search_arima(x,N,d,D,p_max,q_max,P_max,Q_max,Order_max,istationary,s,ic,iapprox,xreg,r,offset,idrift,imean,amethod);
		free(x);
		free(xx);
		free(dx);
		free(diffxreg);
		free(diffdxreg);
		free(diffdx);
		return fit;
	}

	if (N < 10) {
		p_start = (p_start < 1) ? p_start : 1;
		q_start = (q_start < 1) ? q_start : 1;
		P_start = 0;
		Q_start = 0;
	}

	p = p_start = (p_start < p_max) ? p_start : p_max;
	q = q_start = (q_start < q_max) ? q_start : q_max;
	P = P_start = (P_start < P_max) ? P_start : P_max;
	Q = Q_start = (Q_start < Q_max) ? Q_start : Q_max;

	results = (double*)calloc(8*models,sizeof(double));

	bestorder[0] = p;
	bestorder[1] = d;
	bestorder[2] = q;
	bestseasonalorder[0] = P;
	bestseasonalorder[1] = D;
	bestseasonalorder[2] = Q;
	bestseasonalorder[3] = s;

	trace = 0; // Make it variable
	bestconstant = constant;

	bestfit = myarima(x,N,bestorder,bestseasonalorder,bestconstant, ic, trace, iapprox, offset,xreg, r, &amethod);

	best_ic = bestfit->ic;

	set_results(results,1,p,d,q,P,D,Q,constant,bestfit->ic);

	fit->Arima = NULL;
	fit->otype = 1;

	// Null Model with Possible Constant

	order[0] = 0;
	order[1] = d;
	order[2] = 0;
	seasonalorder[0] = 0;
	seasonalorder[1] = D;
	seasonalorder[2] = 0;
	seasonalorder[3] = s;

	fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

	set_results(results,2,0,d,0,0,D,0,constant,fit->myarima->ic);

	if (fit->myarima->ic < bestfit->ic) {
		best_ic = fit->myarima->ic;
		memcpy(bestorder,order,sizeof(int)*3);
		memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
		bestconstant = constant;
		p = q = P = Q = 0;
	}

	k = 2;

	myarima_free(fit->myarima);

	// Basic AR Model

	if (p_max > 0 || P_max > 0) {
		order[0] = (p_max > 0);
		order[1] = d;
		order[2] = 0;
		seasonalorder[0] = (P_max > 0) && (m > 1);
		seasonalorder[1] = D;
		seasonalorder[2] = 0;
		seasonalorder[3] = s;

		fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

		set_results(results,k+1,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

		if (fit->myarima->ic < best_ic) {
			best_ic = fit->myarima->ic;
			memcpy(bestorder,order,sizeof(int)*3);
			memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
			bestconstant = constant;
			p = order[0];
			P = seasonalorder[0];
			q = Q = 0;
		}
		k = k + 1;

		myarima_free(fit->myarima);

	}

	//Basic MA Model

	if ( q_max > 0 || Q_max > 0) {
		order[0] = 0;
		order[1] = d;
		order[2] = (q_max > 0);
		seasonalorder[0] = 0;
		seasonalorder[1] = D;
		seasonalorder[2] = (Q_max > 0) && (m > 1);
		seasonalorder[3] = s;

		fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

		set_results(results,k+1,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

		if (fit->myarima->ic < best_ic) {
			best_ic = fit->myarima->ic;
			memcpy(bestorder,order,sizeof(int)*3);
			memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
			bestconstant = constant;
			q = order[2];
			Q = seasonalorder[2];
			p = P = 0;
		}
		k = k + 1;

		myarima_free(fit->myarima);
	}

	// Null model with no constant

	if (constant) {
		order[0] = 0;
		order[1] = d;
		order[2] = 0;
		seasonalorder[0] = 0;
		seasonalorder[1] = D;
		seasonalorder[2] = 0;
		seasonalorder[3] = s;

		fit->myarima = myarima(x,N,order,seasonalorder,0, ic, trace, iapprox, offset,xreg, r, &amethod);

		set_results(results,k+1,0,d,0,0,D,0,0,fit->myarima->ic);

		if (fit->myarima->ic < best_ic) {
			best_ic = fit->myarima->ic;
			memcpy(bestorder,order,sizeof(int)*3);
			memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
			bestconstant = 0;
			p = q = P = Q = 0;
		}

		k = k + 1;

		myarima_free(fit->myarima);
	}

	startk = 0;

	while (startk < k && k < models) {
		startk = k;

		if (P > 0 && newmodel(p,d,q,P-1,D,Q,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p;
			order[1] = d;
			order[2] = q;
			seasonalorder[0] = P-1;
			seasonalorder[1] = D;
			seasonalorder[2] = Q;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				P = P - 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}

		if (Q > 0 && newmodel(p,d,q,P,D,Q-1,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p;
			order[1] = d;
			order[2] = q;
			seasonalorder[0] = P;
			seasonalorder[1] = D;
			seasonalorder[2] = Q-1;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				Q = Q - 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}

		if (P < P_max && newmodel(p,d,q,P+1,D,Q,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p;
			order[1] = d;
			order[2] = q;
			seasonalorder[0] = P+1;
			seasonalorder[1] = D;
			seasonalorder[2] = Q;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				P = P + 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}

		if (Q < Q_max && newmodel(p,d,q,P,D,Q+1,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p;
			order[1] = d;
			order[2] = q;
			seasonalorder[0] = P;
			seasonalorder[1] = D;
			seasonalorder[2] = Q+1;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				Q = Q + 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}

		if ((Q > 0 ) && (P > 0 )&& newmodel(p,d,q,P-1,D,Q-1,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p;
			order[1] = d;
			order[2] = q;
			seasonalorder[0] = P-1;
			seasonalorder[1] = D;
			seasonalorder[2] = Q-1;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				P = P - 1;
				Q = Q - 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}	

		if ((Q < Q_max ) && (P > 0 )&& newmodel(p,d,q,P-1,D,Q+1,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p;
			order[1] = d;
			order[2] = q;
			seasonalorder[0] = P-1;
			seasonalorder[1] = D;
			seasonalorder[2] = Q+1;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				P = P - 1;
				Q = Q + 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}	

		if ((Q > 0 ) && (P < P_max )&& newmodel(p,d,q,P+1,D,Q-1,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p;
			order[1] = d;
			order[2] = q;
			seasonalorder[0] = P+1;
			seasonalorder[1] = D;
			seasonalorder[2] = Q-1;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				P = P + 1;
				Q = Q - 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}	

		if ((Q < Q_max ) && (P < P_max )&& newmodel(p,d,q,P+1,D,Q+1,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p;
			order[1] = d;
			order[2] = q;
			seasonalorder[0] = P+1;
			seasonalorder[1] = D;
			seasonalorder[2] = Q+1;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				P = P + 1;
				Q = Q + 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}

		if (p > 0 && newmodel(p-1,d,q,P,D,Q,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p-1;
			order[1] = d;
			order[2] = q;
			seasonalorder[0] = P;
			seasonalorder[1] = D;
			seasonalorder[2] = Q;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				p = p - 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}

		if (q > 0 && newmodel(p,d,q-1,P,D,Q,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p;
			order[1] = d;
			order[2] = q-1;
			seasonalorder[0] = P;
			seasonalorder[1] = D;
			seasonalorder[2] = Q;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				q = q - 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}

		if (p < p_max && newmodel(p+1,d,q,P,D,Q,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p+1;
			order[1] = d;
			order[2] = q;
			seasonalorder[0] = P;
			seasonalorder[1] = D;
			seasonalorder[2] = Q;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				p = p + 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}

		if (q < q_max && newmodel(p,d,q+1,P,D,Q,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p;
			order[1] = d;
			order[2] = q+1;
			seasonalorder[0] = P;
			seasonalorder[1] = D;
			seasonalorder[2] = Q;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				q = q + 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}

		if ((q > 0 ) && (p > 0 )&& newmodel(p-1,d,q-1,P,D,Q,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p-1;
			order[1] = d;
			order[2] = q-1;
			seasonalorder[0] = P;
			seasonalorder[1] = D;
			seasonalorder[2] = Q;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				p = p - 1;
				q = q - 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}	

		if ((q < q_max ) && (p > 0 )&& newmodel(p-1,d,q+1,P,D,Q,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p-1;
			order[1] = d;
			order[2] = q+1;
			seasonalorder[0] = P;
			seasonalorder[1] = D;
			seasonalorder[2] = Q;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				p = p - 1;
				q = q + 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}	

		if ((q > 0 ) && (p < p_max )&& newmodel(p+1,d,q-1,P,D,Q,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p+1;
			order[1] = d;
			order[2] = q-1;
			seasonalorder[0] = P;
			seasonalorder[1] = D;
			seasonalorder[2] = Q;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				p = p + 1;
				q = q - 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}	

		if ((q < q_max ) && (p < p_max )&& newmodel(p+1,d,q+1,P,D,Q,constant,results,k)) {
			k = k + 1;
			if (k > models) continue;
			order[0] = p+1;
			order[1] = d;
			order[2] = q+1;
			seasonalorder[0] = P;
			seasonalorder[1] = D;
			seasonalorder[2] = Q;
			seasonalorder[3] = s;

			fit->myarima = myarima(x,N,order,seasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

			set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

			if (fit->myarima->ic < best_ic) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				p = p + 1;
				q = q + 1;
				myarima_free(fit->myarima);
				continue;
			}
			
			myarima_free(fit->myarima);

		}	

		if (idrift || imean) {
			ntconstant = (constant == 1) ? 0 : 1;
			if (newmodel(p,d,q,P,D,Q,ntconstant,results,k)) {
				k = k + 1;
				if (k > models) continue;
				order[0] = p;
				order[1] = d;
				order[2] = q;
				seasonalorder[0] = P;
				seasonalorder[1] = D;
				seasonalorder[2] = Q;
				seasonalorder[3] = s;

				fit->myarima = myarima(x,N,order,seasonalorder,ntconstant, ic, trace, iapprox, offset,xreg, r, &amethod);

				set_results(results,k,order[0],order[1],order[2],seasonalorder[0],seasonalorder[1],seasonalorder[2],constant,fit->myarima->ic);

				if (fit->myarima->ic < best_ic) {
					best_ic = fit->myarima->ic;
					memcpy(bestorder,order,sizeof(int)*3);
					memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
					bestconstant = constant =  ntconstant;
					myarima_free(fit->myarima);
					continue;
				}
				
				myarima_free(fit->myarima);

			}	
		}
	
	}

	printf("DEBUG 5 \n");

	if (k > models) {
		printf("Warning : Stepwise search was stopped early due to reaching the model number limit: %d \n",models);
		k--;
	}

	mdisplay(results,models,8);


	// Delete the previous best model

	myarima_free(bestfit);

	printf("DEBUG 6 k %d \n",k);

	// Refit if Approximation was used

	if (iapprox) {
		icvector = (double*)malloc(sizeof(double)*k);
		icindex = (int*)calloc(k,sizeof(int));

		for(i = 0; i < k;++i) {
			icvector[i] = results[i*8+7];
		}

		sort1d_ascending(icvector,k,icindex);

		for(i = 0; i < k;++i) {
			if (i > 0) {
				myarima_free(fit->myarima);
			}
			
			printf(" Free %d \n",i+1);
			order[0] = results[icindex[i]*8];
			order[1] = d;
			order[2] = results[icindex[i]*8 + 2];
			seasonalorder[0] = results[icindex[i]*8 + 3];
			seasonalorder[1] = D;
			seasonalorder[2] = results[icindex[i]*8 + 5];
			seasonalorder[3] = s;
			fit->myarima = myarima(x,N,order,seasonalorder,results[icindex[i]*8 + 6], ic, trace, 0, offset,xreg, r, &amethod);

			if (fit->myarima->ic < DBL_MAX) {
				best_ic = fit->myarima->ic;
				memcpy(bestorder,order,sizeof(int)*3);
				memcpy(bestseasonalorder,seasonalorder,sizeof(int)*4);
				bestconstant = constant;
				free(icvector);
				free(icindex);
				iapprox = 0;
				break;
			}
		}

		if (iapprox != 0) {
			free(icvector);
			free(icindex);
		}
		
	}

	// Refit The Best Model

	fit->myarima = myarima(x,N,bestorder,bestseasonalorder,constant, ic, trace, iapprox, offset,xreg, r, &amethod);

	

	free(x);
	free(xx);
	free(dx);
	free(diffxreg);
	free(diffdxreg);
	free(diffdx);
	free(results);
	return fit;
}

void ar_exec(ar_object obj, double *inp) {
	int p, N;
	double loglik;
	double *hess;
	N = obj->N;

	obj->p = ar_estimate(inp, N, obj->method);
	//printf("\n p %d \n", obj->p);
	p = obj->p;
	int cssml = 0;

	if (obj->method == 2) {
		hess = (double*)malloc(sizeof(double)* (p + 1) * (p + 1));
		obj->retval = as154(inp, obj->N, obj->optmethod, obj->p, 0, 0 , obj->params, NULL, &obj->mean, &obj->var, obj->params + p, &loglik,
			hess,cssml);
		loglik = -0.5 * (N * (2 * loglik + 1.0 + log(2 * 3.14159)));
		obj->aic = -2.0 * loglik + 2.0 * (obj->p + 1) + 2.0;
		free(hess);
	}
	else if (obj->method == 0) {
		obj->mean = mean(inp, N);
		ywalg2(inp, N, p, obj->params, &obj->var);
	}
	else if (obj->method == 1) {
		obj->mean = mean(inp, N);
		burgalg(inp, N - 1, p, obj->params, &obj->var);
	}

	obj->order = obj->p;


}

void arima_predict(arima_object obj,double *inp, int L,double *xpred,double *amse) {
	int d,i,N,ip,iq,ir;
	double *delta,*W,*resid,*phi,*theta;
	double wmean;

	d = obj->d;
	N = obj->N;
	ip = obj->p;
	iq = obj->q;
	delta = (double*)malloc(sizeof(double)* (d + 1));
	W = (double*)malloc(sizeof(double)* N);
	resid = (double*)malloc(sizeof(double)* N);

	ir = iq + 1;

	if (ir < ip) {
		ir = ip;
	}

	phi = (double*)malloc(sizeof(double)* ir);
	theta = (double*)malloc(sizeof(double)* ir);
	wmean = 0.0;
	if (d > 0) {
		deld(d, delta);
	}
	else {
		*delta = 1.0;
		wmean = obj->mean;
	}
	for (i = 1; i <= d; ++i) {
		delta[i] = -1.0 * delta[i];
	}
	for (i = 0; i < N; ++i) {
		W[i] = inp[i];
		if (d == 0) {
			W[i] -= wmean;
		}
		resid[i] = obj->res[i];
	}
	for (i = 0; i < ir; ++i) {
		phi[i] = theta[i] = 0.0;
	}

	for (i = 0; i < ip; ++i) {
		phi[i] = obj->phi[i];
	}

	for (i = 0; i < iq; ++i) {
		theta[i] = -1.0 *  obj->theta[i];
	}

	forkal(ip, iq, d, phi, theta, delta+1, N, W, resid, L, xpred, amse);

	for (i = 0; i < L; ++i) {
		xpred[i] += wmean;
	}

	free(delta);
	free(W);
	free(resid);
	free(phi);
	free(theta);
}

void ar_predict(ar_object obj, double *inp, int L, double *xpred, double *amse) {
	int d, i, N, ip, iq, ir;
	double *delta, *W, *resid, *phi, *theta;
	double wmean;

	d = 0;
	N = obj->N;
	ip = obj->p;
	iq = 0;
	delta = (double*)malloc(sizeof(double)* (d + 1));
	W = (double*)malloc(sizeof(double)* N);
	resid = (double*)malloc(sizeof(double)* N);

	ir = iq + 1;

	if (ir < ip) {
		ir = ip;
	}

	phi = (double*)malloc(sizeof(double)* ir);
	theta = (double*)malloc(sizeof(double)* ir);
	wmean = 0.0;
	if (d > 0) {
		deld(d, delta);
	}
	else {
		*delta = 1.0;
		wmean = obj->mean;
	}
	for (i = 1; i <= d; ++i) {
		delta[i] = -1.0 * delta[i];
	}
	for (i = 0; i < N; ++i) {
		W[i] = inp[i];
		if (d == 0) {
			W[i] -= wmean;
		}
		resid[i] = obj->res[i];
	}
	for (i = 0; i < ir; ++i) {
		phi[i] = theta[i] = 0.0;
	}

	for (i = 0; i < ip; ++i) {
		phi[i] = obj->phi[i];
	}

	for (i = 0; i < iq; ++i) {
		theta[i] = -1.0 *  0.0;
	}

	forkal(ip, iq, d, phi, theta, delta + 1, N, W, resid, L, xpred, amse);

	for (i = 0; i < L; ++i) {
		xpred[i] += wmean;
	}

	free(delta);
	free(W);
	free(resid);
	free(phi);
	free(theta);
}

void ar(double *inp, int N, int p,int method,double *phi,double *var) {
	/*
	Method 0 : Yule Walker
	Method 1 : Burg
	Method 2 : MLE
	*/

	double wmean,loglik;
	double *resid,*hess;
	int cssml = 0;

	if (method == 0) {
		ywalg2(inp, N, p, phi,var);
	}
	else if (method == 1) {
		burgalg(inp, N - 1, p, phi, var);
	}
	else if (method == 2) {
		resid = (double*)malloc(sizeof(double)* N);
		hess = (double*)malloc(sizeof(double)* (p+1) * (p+1));
		as154(inp,N,7,p,0,0,phi,NULL, &wmean, var, resid, &loglik,hess,cssml);
		free(resid);
		free(hess);
	}
	else {
		printf("The program only accepts 0,1 and 2 as input values for method");
	}

}

void arima_setMethod(arima_object obj, int value) {
	if (value == 0) {
		obj->method = 0;
	}
	else if (value == 1) {
		obj->method = 1;
	}
	else if (value == 2) {
		obj->method = 2;
	}
	else {
		printf("\n Acceptable Numerical Values 0 - MLE, 1 - CSS, 2 - Box-Jenkins \n");
	}
}

void sarima_setMethod(sarima_object obj, int value) {
	if (value == 0) {
		obj->method = 0;
	}
	else if (value == 1) {
		obj->method = 1;
	}
	else if (value == 2) {
		obj->method = 2;
	}
	else {
		printf("\n Acceptable Numerical Values 0 - MLE, 1 - CSS, 2 - Box-Jenkins \n");
	}
}

void sarimax_setMethod(sarimax_object obj, int value) {
	if (value == 0) {
		obj->method = 0;
	}
	else if (value == 1) {
		obj->method = 1;
	}
	else if (value == 2) {
		obj->method = 2;
	}
	else {
		printf("\n Acceptable Numerical Values 0 - CSS-MLE, 1 - MLE, 2 - CSS \n");
	}
}

void auto_arima_setMethod(auto_arima_object obj, int value) {
	if (value == 0) {
		obj->method = 0;
	}
	else if (value == 1) {
		obj->method = 1;
	}
	else if (value == 2) {
		obj->method = 2;
	}
	else {
		printf("\n Acceptable Numerical Values 0 - CSS-MLE, 1 - MLE, 2 - CSS \n");
	}
}


void sarimax_setParams(sarimax_object obj, double *phi, double *theta, double *PHI, double *THETA) {
	int i;
	if (phi) {
		for (i = 0; i < obj->p;++i) {
			obj->phi[i] = phi[i];
		}
	} else {
		for (i = 0; i < obj->p;++i) {
			obj->phi[i] = 0.0;
		}
	}

	if (theta) {
		for (i = 0; i < obj->q;++i) {
			obj->theta[i] = theta[i];
		}
	} else {
		for (i = 0; i < obj->q;++i) {
			obj->theta[i] = 0.0;
		}
	}

	if (PHI) {
		for (i = 0; i < obj->P;++i) {
			obj->PHI[i] = PHI[i];
		}
	} else {
		for (i = 0; i < obj->P;++i) {
			obj->PHI[i] = 0.0;
		}
	}

	if (THETA) {
		for (i = 0; i < obj->Q;++i) {
			obj->THETA[i] = THETA[i];
		}
	} else {
		for (i = 0; i < obj->Q;++i) {
			obj->THETA[i] = 0.0;
		}
	}

	obj->start = 1;
}

void arima_setCSSML(arima_object obj, int cssml) {
	/*
	Uses CSS before MLE if cssml = 1.
	Only uses MLE if cssml = 0
	Only applicable with method 0
	*/

	if (cssml == 1) {
		obj->cssml = 1;
	} else if (cssml == 0) {
		obj->cssml = 0;
	} else {
		printf("cssml only accepts two values 1 and 0 \n");
		exit(-1);
	}

}

void sarima_setCSSML(sarima_object obj, int cssml) {
	/*
	Uses CSS before MLE if cssml = 1.
	Only uses MLE if cssml = 0
	Only applicable with method 0
	*/

	if (cssml == 1) {
		obj->cssml = 1;
	} else if (cssml == 0) {
		obj->cssml = 0;
	} else {
		printf("cssml only accepts two values 1 and 0 \n");
		exit(-1);
	}

}

void arima_setOptMethod(arima_object obj, int value) {
	/*
	* Method 0 - Nelder-Mead
	* Method 1 - Newton Line Search
	* Method 2 - Newton Trust Region - Hook Step
	* Method 3 - Newton Trust Region - Double Dog-Leg
	* Method 4 - Conjugate Gradient
	* Method 5 - BFGS
	* Method 6 - Limited Memory BFGS
	* Method 7 - BFGS More-Thuente Line Search 
	*/
	if (value == 0) {
		obj->optmethod = 0;
	}
	else if (value == 1) {
		obj->optmethod = 1;
	}
	else if (value == 2) {
		obj->optmethod = 2;
	}
	else if (value == 3) {
		obj->optmethod = 3;
	}
	else if (value == 4) {
		obj->optmethod = 4;
	}
	else if (value == 5) {
		obj->optmethod = 5;
	}
	else if (value == 6) {
		obj->optmethod = 6;
	}
	else if (value == 7) {
		obj->optmethod = 7;
	}
	else {
		printf("\n Acceptable Numerical Values 0,1,2,3,4,5,6,7 \n");
		printf("\n Method 0 - Nelder-Mead \n");
		printf("\n Method 1 - Newton Line Search \n");
		printf("\n Method 2 - Newton Trust Region - Hook Step \n");
		printf("\n Method 3 - Newton Trust Region - Double Dog-Leg \n");
		printf("\n Method 4 - Conjugate Gradient \n");
		printf("\n Method 5 - BFGS \n");
		printf("\n Method 6 - Limited Memory BFGS \n");
		printf("\n Method 7 - BFGS More-Thuente Line Search \n");
	}
}

void sarima_setOptMethod(sarima_object obj, int value) {
	/*
	* Method 0 - Nelder-Mead
	* Method 1 - Newton Line Search
	* Method 2 - Newton Trust Region - Hook Step
	* Method 3 - Newton Trust Region - Double Dog-Leg
	* Method 4 - Conjugate Gradient
	* Method 5 - BFGS
	* Method 6 - Limited Memory BFGS
	* Method 7 - BFGS More-Thuente Line Search
	*/
	if (value == 0) {
		obj->optmethod = 0;
	}
	else if (value == 1) {
		obj->optmethod = 1;
	}
	else if (value == 2) {
		obj->optmethod = 2;
	}
	else if (value == 3) {
		obj->optmethod = 3;
	}
	else if (value == 4) {
		obj->optmethod = 4;
	}
	else if (value == 5) {
		obj->optmethod = 5;
	}
	else if (value == 6) {
		obj->optmethod = 6;
	}
	else if (value == 7) {
		obj->optmethod = 7;
	}
	else {
		printf("\n Acceptable Numerical Values 0,1,2,3,4,5,6,7 \n");
		printf("\n Method 0 - Nelder-Mead \n");
		printf("\n Method 1 - Newton Line Search \n");
		printf("\n Method 2 - Newton Trust Region - Hook Step \n");
		printf("\n Method 3 - Newton Trust Region - Double Dog-Leg \n");
		printf("\n Method 4 - Conjugate Gradient \n");
		printf("\n Method 5 - BFGS \n");
		printf("\n Method 6 - Limited Memory BFGS \n");
		printf("\n Method 7 - BFGS More-Thuente Line Search \n");
	}
}

void sarimax_setOptMethod(sarimax_object obj, int value) {
	/*
	* Method 0 - Nelder-Mead
	* Method 1 - Newton Line Search
	* Method 2 - Newton Trust Region - Hook Step
	* Method 3 - Newton Trust Region - Double Dog-Leg
	* Method 4 - Conjugate Gradient
	* Method 5 - BFGS
	* Method 6 - Limited Memory BFGS
	* Method 7 - BFGS More-Thuente Line Search
	*/
	if (value == 0) {
		obj->optmethod = 0;
	}
	else if (value == 1) {
		obj->optmethod = 1;
	}
	else if (value == 2) {
		obj->optmethod = 2;
	}
	else if (value == 3) {
		obj->optmethod = 3;
	}
	else if (value == 4) {
		obj->optmethod = 4;
	}
	else if (value == 5) {
		obj->optmethod = 5;
	}
	else if (value == 6) {
		obj->optmethod = 6;
	}
	else if (value == 7) {
		obj->optmethod = 7;
	}
	else {
		printf("\n Acceptable Numerical Values 0,1,2,3,4,5,6,7 \n");
		printf("\n Method 0 - Nelder-Mead \n");
		printf("\n Method 1 - Newton Line Search \n");
		printf("\n Method 2 - Newton Trust Region - Hook Step \n");
		printf("\n Method 3 - Newton Trust Region - Double Dog-Leg \n");
		printf("\n Method 4 - Conjugate Gradient \n");
		printf("\n Method 5 - BFGS \n");
		printf("\n Method 6 - Limited Memory BFGS \n");
		printf("\n Method 7 - BFGS More-Thuente Line Search \n");
	}
}

void auto_arima_setOptMethod(auto_arima_object obj, int value) {
	/*
	* Method 0 - Nelder-Mead
	* Method 1 - Newton Line Search
	* Method 2 - Newton Trust Region - Hook Step
	* Method 3 - Newton Trust Region - Double Dog-Leg
	* Method 4 - Conjugate Gradient
	* Method 5 - BFGS
	* Method 6 - Limited Memory BFGS
	* Method 7 - BFGS More-Thuente Line Search
	*/
	if (value == 0) {
		obj->optmethod = 0;
	}
	else if (value == 1) {
		obj->optmethod = 1;
	}
	else if (value == 2) {
		obj->optmethod = 2;
	}
	else if (value == 3) {
		obj->optmethod = 3;
	}
	else if (value == 4) {
		obj->optmethod = 4;
	}
	else if (value == 5) {
		obj->optmethod = 5;
	}
	else if (value == 6) {
		obj->optmethod = 6;
	}
	else if (value == 7) {
		obj->optmethod = 7;
	}
	else {
		printf("\n Acceptable Numerical Values 0,1,2,3,4,5,6,7 \n");
		printf("\n Method 0 - Nelder-Mead \n");
		printf("\n Method 1 - Newton Line Search \n");
		printf("\n Method 2 - Newton Trust Region - Hook Step \n");
		printf("\n Method 3 - Newton Trust Region - Double Dog-Leg \n");
		printf("\n Method 4 - Conjugate Gradient \n");
		printf("\n Method 5 - BFGS \n");
		printf("\n Method 6 - Limited Memory BFGS \n");
		printf("\n Method 7 - BFGS More-Thuente Line Search \n");
	}
}

void arima_vcov(arima_object obj, double *vcov) {
	int i;

	for (i = 0; i < obj->lvcov; ++i) {
		vcov[i] = obj->vcov[i];
	}
}

void sarima_vcov(sarima_object obj, double *vcov) {
	int i;

	for (i = 0; i < obj->lvcov; ++i) {
		vcov[i] = obj->vcov[i];
	}
}

void sarimax_vcov(sarimax_object obj, double *vcov) {
	int i;

	for (i = 0; i < obj->lvcov; ++i) {
		vcov[i] = obj->vcov[i];
	}
}

void auto_arima_setApproximation(auto_arima_object obj, int approximation) {
	if (approximation == 0 || approximation == 1) {
		obj->approximation = approximation;
	} else {
		printf("Approximation parameter accepts only two values - 0 or 1 \n");
	}
}

void auto_arima_setStepwise(auto_arima_object obj, int stepwise) {
	if (stepwise == 0 || stepwise == 1) {
		obj->stepwise = stepwise;
	} else {
		printf("Stepwise parameter accepts only two values - 0 or 1 \n");
	}
}

void auto_arima_setStationary(auto_arima_object obj, int stationary) {
	if (stationary == 0 || stationary == 1) {
		obj->stationary = stationary;
	} else {
		printf("Stationary parameter accepts only two values - 0 or 1 \n");
	}
}

void auto_arima_setSeasonal(auto_arima_object obj, int seasonal) {
	if (seasonal == 0 || seasonal == 1) {
		obj->seasonal = seasonal;
	} else {
		printf("Seasonal parameter accepts only two values - 0 or 1 \n");
	}
}

void auto_arima_setStationarityParameters(auto_arima_object obj,const char *test, double alpha, const char *type) {
	if (!strcmp(test,"kpss") || !strcmp(test,"pp") || !strcmp(test,"df") || !strcmp(test,"adf")) {
		strcpy(obj->test,test);
	} else {
		printf("Only three tests are allowed - kpss, df and pp \n");
		exit(-1);
	}

	if (!strcmp(type,"level") || !strcmp(type,"trend")) {
		strcpy(obj->type,type);
	} else {
		printf("Only two tests are allowed - level and trend \n");
		exit(-1);
	}

	obj->alpha_test = alpha;
}

void auto_arima_setSeasonalParameters(auto_arima_object obj,const char *test, double alpha) {
	if (!strcmp(test,"ocsb") || !strcmp(test,"seas")) {
		strcpy(obj->seas,test);
	} else {
		printf("Only two tests are allowed - seas and ocsb \n");
		exit(-1);
	}

	obj->alpha_seas = alpha;
}

void auto_arima_allowMean(auto_arima_object obj, int allowmean) {
	if (allowmean == 0 || allowmean == 1) {
		obj->imean = allowmean;
	} else {
		printf("Allowmean parameter accepts only two values - 0 or 1 \n");
	}
}

void auto_arima_allowDrift(auto_arima_object obj, int allowdrift) {
	if (allowdrift == 0 || allowdrift == 1) {
		obj->idrift = allowdrift;
	} else {
		printf("Allowdrift parameter accepts only two values - 0 or 1 \n");
	}
}

void sarima_predict(sarima_object obj, double *inp, int L, double *xpred, double *amse) {
	int d, i, N, ip, iq, ir,D,P,Q,s,p,q,t,ps,qs,j;
	double *coef1,*coef2,*delta, *W, *resid, *phi, *theta;
	double wmean;

	d = obj->d;
	N = obj->N;
	p = obj->p;
	q = obj->q;
	D = obj->D;
	P = obj->P;
	Q = obj->Q;
	s = obj->s;

	ip = p + s * P;
	iq = q + s * Q;
	ir = p + s * P;

	t = 1 + q + s*Q;
	if (ir < t) {
		ir = t;
	}
	ps = P;
	qs = Q;
	coef1 = (double*)malloc(sizeof(double)* (d + 1));
	coef2 = (double*)malloc(sizeof(double)* (D*s + 1));
	delta = (double*)malloc(sizeof(double)* (d + D*s + 1));
	W = (double*)malloc(sizeof(double)* N);
	resid = (double*)malloc(sizeof(double)* N);


	phi = (double*)malloc(sizeof(double)* ir);
	theta = (double*)malloc(sizeof(double)* ir);
	wmean = 0.0;
	coef1[0] = coef2[0] = 1.0;

	if (d == 0 && D == 0) {
		*delta = 1.0;
		wmean = obj->mean;
	}

	if (d > 0) {
		deld(d, coef1);
	}

	if (D > 0) {
		delds(D, s, coef2);
	}

	conv(coef1, d + 1, coef2, D*s + 1, delta);
	
	for (i = 1; i <= d+D*s; ++i) {
		delta[i] = -1.0 * delta[i];
	}
	for (i = 0; i < N; ++i) {
		W[i] = inp[i];
		if (d == 0 && D == 0) {
			W[i] -= wmean;
		}
		resid[i] = obj->res[i];
	}
	for (i = 0; i < ir; ++i) {
		phi[i] = theta[i] = 0.0;
	}

	for (i = 0; i < p; ++i) {
		phi[i] = obj->phi[i];
	}

	for (i = 0; i < q; ++i) {
		theta[i] = -1.0 *  obj->theta[i];
	}

	for (j = 0; j < ps; ++j) {
		phi[(j + 1)*s - 1] += obj->PHI[j];
		for (i = 0; i < p; ++i) {
			phi[(j + 1)*s + i] -= obj->phi[i] * obj->PHI[j];
		}
	}

	for (j = 0; j < qs; ++j) {
		theta[(j + 1)*s - 1] -= obj->THETA[j];
		for (i = 0; i < q; ++i) {
			theta[(j + 1)*s + i] += obj->theta[i] * obj->THETA[j];
		}
	}

	forkal(ip, iq, d+D*s, phi, theta, delta + 1, N, W, resid, L, xpred, amse);

	for (i = 0; i < L; ++i) {
		xpred[i] += wmean;
	}

	free(coef1);
	free(coef2);
	free(delta);
	free(W);
	free(resid);
	free(phi);
	free(theta);
}

void sarimax_predict(sarimax_object obj, double *inp, double *xreg, int L,double *newxreg, double *xpred, double *amse) {
	int d, i, N, ip, iq, ir,D,P,Q,s,p,q,t,ps,qs,j;
	double *coef1,*coef2,*delta, *W, *resid, *phi, *theta;
	double wmean;

	d = obj->d;
	N = obj->N;
	p = obj->p;
	q = obj->q;
	D = obj->D;
	P = obj->P;
	Q = obj->Q;
	s = obj->s;

	ip = p + s * P;
	iq = q + s * Q;
	ir = p + s * P;

	t = 1 + q + s*Q;
	if (ir < t) {
		ir = t;
	}
	ps = P;
	qs = Q;
	coef1 = (double*)malloc(sizeof(double)* (d + 1));
	coef2 = (double*)malloc(sizeof(double)* (D*s + 1));
	delta = (double*)malloc(sizeof(double)* (d + D*s + 1));
	W = (double*)malloc(sizeof(double)* N);
	resid = (double*)malloc(sizeof(double)* N);


	phi = (double*)malloc(sizeof(double)* ir);
	theta = (double*)malloc(sizeof(double)* ir);
	wmean = 0.0;
	coef1[0] = coef2[0] = 1.0;

	if (d == 0 && D == 0) {
		*delta = 1.0;
		wmean = obj->mean;
	}

	if (d > 0) {
		deld(d, coef1);
	}

	if (D > 0) {
		delds(D, s, coef2);
	}

	conv(coef1, d + 1, coef2, D*s + 1, delta);
	
	for (i = 1; i <= d+D*s; ++i) {
		delta[i] = -1.0 * delta[i];
	}
	for (i = 0; i < N; ++i) {
		W[i] = inp[i];
		if (d == 0 && D == 0) {
			W[i] -= wmean;
		}
		if (obj->r > 0) {
			for(j = 0; j < obj->r;++j) {
				W[i] -= obj->exog[j] * xreg[j*N+i];
			}
		}
		resid[i] = obj->res[i];
	}
	for (i = 0; i < ir; ++i) {
		phi[i] = theta[i] = 0.0;
	}

	for (i = 0; i < p; ++i) {
		phi[i] = obj->phi[i];
	}

	for (i = 0; i < q; ++i) {
		theta[i] = -1.0 *  obj->theta[i];
	}

	for (j = 0; j < ps; ++j) {
		phi[(j + 1)*s - 1] += obj->PHI[j];
		for (i = 0; i < p; ++i) {
			phi[(j + 1)*s + i] -= obj->phi[i] * obj->PHI[j];
		}
	}

	for (j = 0; j < qs; ++j) {
		theta[(j + 1)*s - 1] -= obj->THETA[j];
		for (i = 0; i < q; ++i) {
			theta[(j + 1)*s + i] += obj->theta[i] * obj->THETA[j];
		}
	}

	forkal(ip, iq, d+D*s, phi, theta, delta + 1, N, W, resid, L, xpred, amse);

	for (i = 0; i < L; ++i) {
		xpred[i] += wmean;
		for(j = 0; j < obj->r;++j) {
			xpred[i] += obj->exog[j] * newxreg[j*L+i];
		}
	}

	free(coef1);
	free(coef2);
	free(delta);
	free(W);
	free(resid);
	free(phi);
	free(theta);
}

void auto_arima_predict(auto_arima_object obj, double *inp, double *xreg, int L,double *newxreg, double *xpred, double *amse) {
	int d, i, N, ip, iq, ir,D,P,Q,s,p,q,t,ps,qs,j,diter;
	double *coef1,*coef2,*delta, *W, *resid, *phi, *theta;
	double wmean;

	d = obj->d;
	N = obj->N;
	p = obj->p;
	q = obj->q;
	D = obj->D;
	P = obj->P;
	Q = obj->Q;
	s = obj->s;

	ip = p + s * P;
	iq = q + s * Q;
	ir = p + s * P;

	t = 1 + q + s*Q;
	if (ir < t) {
		ir = t;
	}
	ps = P;
	qs = Q;
	coef1 = (double*)malloc(sizeof(double)* (d + 1));
	coef2 = (double*)malloc(sizeof(double)* (D*s + 1));
	delta = (double*)malloc(sizeof(double)* (d + D*s + 1));
	W = (double*)malloc(sizeof(double)* N);
	resid = (double*)malloc(sizeof(double)* N);

	diter = 0;

	if (obj->idrift == 1) {
		diter = 1;		
	}

	phi = (double*)malloc(sizeof(double)* ir);
	theta = (double*)malloc(sizeof(double)* ir);
	wmean = 0.0;
	coef1[0] = coef2[0] = 1.0;

	if (d == 0 && D == 0) {
		*delta = 1.0;
		wmean = obj->mean;
	}

	if (d > 0) {
		deld(d, coef1);
	}

	if (D > 0) {
		delds(D, s, coef2);
	}

	conv(coef1, d + 1, coef2, D*s + 1, delta);
	
	for (i = 1; i <= d+D*s; ++i) {
		delta[i] = -1.0 * delta[i];
	}
	for (i = 0; i < N; ++i) {
		W[i] = inp[i];
		if (d == 0 && D == 0) {
			W[i] -= wmean;
		}
		if (obj->idrift == 1) {
			W[i] -= obj->exog[0]*(i+1);
		}
		if (obj->sarimax->r > 0) {
			for(j = diter; j < obj->r;++j) {
				W[i] -= obj->exog[j] * xreg[(j-diter)*N+i];
			}
		}
		resid[i] = obj->res[i];
	}
	for (i = 0; i < ir; ++i) {
		phi[i] = theta[i] = 0.0;
	}

	for (i = 0; i < p; ++i) {
		phi[i] = obj->phi[i];
	}

	for (i = 0; i < q; ++i) {
		theta[i] = -1.0 *  obj->theta[i];
	}

	for (j = 0; j < ps; ++j) {
		phi[(j + 1)*s - 1] += obj->PHI[j];
		for (i = 0; i < p; ++i) {
			phi[(j + 1)*s + i] -= obj->phi[i] * obj->PHI[j];
		}
	}

	for (j = 0; j < qs; ++j) {
		theta[(j + 1)*s - 1] -= obj->THETA[j];
		for (i = 0; i < q; ++i) {
			theta[(j + 1)*s + i] += obj->theta[i] * obj->THETA[j];
		}
	}

	forkal(ip, iq, d+D*s, phi, theta, delta + 1, N, W, resid, L, xpred, amse);

	for (i = 0; i < L; ++i) {
		xpred[i] += wmean;
		if (obj->idrift == 1) {
			xpred[i] += obj->exog[0]*(i+1+N);
		}
		for(j = diter; j < obj->r;++j) {
			xpred[i] += obj->exog[j] * newxreg[(j - diter)*L+i];
		}
	}

	free(coef1);
	free(coef2);
	free(delta);
	free(W);
	free(resid);
	free(phi);
	free(theta);
}


void sarimax_wrapper_predict(sarimax_wrapper_object obj, double *inp, double *xreg, int L,double *newxreg, double *xpred, double *amse) {
	int d, i, N, ip, iq, ir,D,P,Q,s,p,q,t,ps,qs,j,diter;
	double *coef1,*coef2,*delta, *W, *resid, *phi, *theta;
	double wmean;

	d = obj->sarimax->d;
	N = obj->sarimax->N;
	p = obj->sarimax->p;
	q = obj->sarimax->q;
	D = obj->sarimax->D;
	P = obj->sarimax->P;
	Q = obj->sarimax->Q;
	s = obj->sarimax->s;

	ip = p + s * P;
	iq = q + s * Q;
	ir = p + s * P;

	t = 1 + q + s*Q;
	if (ir < t) {
		ir = t;
	}
	ps = P;
	qs = Q;
	coef1 = (double*)malloc(sizeof(double)* (d + 1));
	coef2 = (double*)malloc(sizeof(double)* (D*s + 1));
	delta = (double*)malloc(sizeof(double)* (d + D*s + 1));
	W = (double*)malloc(sizeof(double)* N);
	resid = (double*)malloc(sizeof(double)* N);

	diter = 0;

	if (obj->idrift == 1) {
		diter = 1;		
	}


	phi = (double*)malloc(sizeof(double)* ir);
	theta = (double*)malloc(sizeof(double)* ir);
	wmean = 0.0;
	coef1[0] = coef2[0] = 1.0;

	if (d == 0 && D == 0) {
		*delta = 1.0;
		wmean = obj->sarimax->mean;
	}

	if (d > 0) {
		deld(d, coef1);
	}

	if (D > 0) {
		delds(D, s, coef2);
	}

	conv(coef1, d + 1, coef2, D*s + 1, delta);
	
	for (i = 1; i <= d+D*s; ++i) {
		delta[i] = -1.0 * delta[i];
	}
	for (i = 0; i < N; ++i) {
		W[i] = inp[i];
		if (d == 0 && D == 0) {
			W[i] -= wmean;
		}
		if (obj->idrift == 1) {
			W[i] -= obj->sarimax->exog[0]*(i+1);
		}
		if (obj->sarimax->r > 0) {
			for(j = diter; j < obj->sarimax->r;++j) {
				W[i] -= obj->sarimax->exog[j] * xreg[(j-diter)*N+i];
			}
		}
		resid[i] = obj->sarimax->res[i];
	}
	for (i = 0; i < ir; ++i) {
		phi[i] = theta[i] = 0.0;
	}

	for (i = 0; i < p; ++i) {
		phi[i] = obj->sarimax->phi[i];
	}

	for (i = 0; i < q; ++i) {
		theta[i] = -1.0 *  obj->sarimax->theta[i];
	}

	for (j = 0; j < ps; ++j) {
		phi[(j + 1)*s - 1] += obj->sarimax->PHI[j];
		for (i = 0; i < p; ++i) {
			phi[(j + 1)*s + i] -= obj->sarimax->phi[i] * obj->sarimax->PHI[j];
		}
	}

	for (j = 0; j < qs; ++j) {
		theta[(j + 1)*s - 1] -= obj->sarimax->THETA[j];
		for (i = 0; i < q; ++i) {
			theta[(j + 1)*s + i] += obj->sarimax->theta[i] * obj->sarimax->THETA[j];
		}
	}

	forkal(ip, iq, d+D*s, phi, theta, delta + 1, N, W, resid, L, xpred, amse);

	for (i = 0; i < L; ++i) {
		xpred[i] += wmean;
		if (obj->idrift == 1) {
			xpred[i] += obj->sarimax->exog[0]*(i+1+N);
		}
		for(j = diter; j < obj->sarimax->r;++j) {
			xpred[i] += obj->sarimax->exog[j] * newxreg[(j - diter)*L+i];
		}
	}

	free(coef1);
	free(coef2);
	free(delta);
	free(W);
	free(resid);
	free(phi);
	free(theta);
}

void sarima_summary(sarima_object obj) {
	int i,pq,t;
	pq = obj->p + obj->q + obj->P + obj->Q + obj->M;
	if (obj->method == 0 || obj->method == 1) {
		printf("\n\n Exit Status \n");
		printf("Return Code : %d \n", obj->retval);
		printf("Exit Message : ");

		if (obj->retval == 0) {
			printf("Input Error");
		}
		else if (obj->retval == 1) {
			printf("Probable Success");
		}
		else if (obj->retval == 4) {
			printf("Optimization Routine didn't converge");
		}
		else if (obj->retval == 15) {
			printf("Optimization Routine Encountered Inf/Nan Values");
		}
	}
	printf("\n\n");
	printf("  ARIMA Seasonal Order : ( %d, %d, %d) * (%d, %d, %d) \n",obj->p,obj->d,obj->q, obj->P,obj->D,obj->Q );
	printf("\n");

	printf("%-20s%-20s%-20s \n\n", "Coefficients", "Value", "Standard Error");
	for (i = 0; i < obj->p; ++i) {
		printf("AR%-15d%-20g%-20g \n", i + 1, obj->phi[i], sqrt(obj->vcov[i + pq*i]));
	}
	for (i = 0; i < obj->q; ++i) {
		t = obj->p + i;
		printf("MA%-15d%-20g%-20g \n", i + 1, obj->theta[i], sqrt(obj->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->P; ++i) {
		t = obj->p + obj->q + i;
		printf("SAR%-14d%-20g%-20g \n", i + 1, obj->PHI[i], sqrt(obj->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->Q; ++i) {
		t = obj->p + obj->q + obj->P + i;
		printf("SMA%-14d%-20g%-20g \n", i + 1, obj->THETA[i], sqrt(obj->vcov[t + pq * t]));
	}
	printf("\n");
	t = obj->p + obj->q + obj->P + obj->Q;
	if (obj->M == 1) {
		printf("%-17s%-20g%-20g \n", "MEAN", obj->mean, sqrt(obj->vcov[t + pq * t]));
	}
	else {
		printf("%-17s%-20g \n", "MEAN", obj->mean);
	}
	printf("\n");
	printf("%-17s%-20g \n", "SIGMA^2", obj->var);
	printf("\n");
	printf("ESTIMATION METHOD : ");
	if (obj->method == 0) {
		printf("MLE");
	}
	else if (obj->method == 1) {
		printf("CSS");
	}
	else if (obj->method == 2) {
		printf("BOX-JENKINS");
	}
	printf("\n\n");
	printf("OPTIMIZATION METHOD : ");
	if (obj->method == 2) {
		printf("Newton-Raphson");
	}
	else {
		if (obj->optmethod == 0) {
			printf("Nelder-Mead");
		}
		else if (obj->optmethod == 1) {
			printf("Newton Line Search");
		}
		else if (obj->optmethod == 2) {
			printf("Newton Trust Region - Hook Step");
		}
		else if (obj->optmethod == 3) {
			printf("Newton Trust Region - Double Dog-Leg");
		}
		else if (obj->optmethod == 4) {
			printf("Conjugate Gradient");
		}
		else if (obj->optmethod == 5) {
			printf("BFGS");
		}
		else if (obj->optmethod == 6) {
			printf("L-BFGS");
		}
		else if (obj->optmethod == 7) {
			printf("BFGS More-Thuente Line Search");
		}
	}
	printf("\n\n");
	if (obj->method == 0 || obj->method == 1) {
		printf("Log Likelihood : %g ", obj->loglik);
		printf("\n\n");
	}
	else {
		printf("Log Likelihood : Unavailable ");
		printf("\n\n");
	}
	if (obj->method == 0) {
		printf("AIC criterion : %g ", obj->aic);
		printf("\n\n");
	}
	else {
		printf("AIC Criterion : Unavailable ");
		printf("\n\n");
	}
	/*
	printf("EQUATION FORM : x[t] ");
	for (i = 0; i < obj->p; ++i) {
		if (obj->phi[i] > 0) {
			printf("- %g*x[t - %d%s", fabs(obj->phi[i]), i + 1, "]");
		}
		else if (obj->phi[i] < 0) {
			printf("+ %g*x[t - %d%s", fabs(obj->phi[i]), i + 1, "] ");
		}
	}

	printf("=");
	if (obj->mean != 0.0) {
		printf(" %g +", obj->mean);
	}
	printf(" a[t] ");
	for (i = 0; i < obj->q; ++i) {
		if (obj->theta[i] > 0) {
			printf("- %g*a[t - %d%s", fabs(obj->theta[i]), i + 1, "]");
		}
		else if (obj->phi[i] < 0) {
			printf("+ %g*a[t - %d%s", fabs(obj->theta[i]), i + 1, "] ");
		}
	}
	printf("\n\n");
	*/
}

void arima_summary(arima_object obj) {
	int i, pq,t;
	pq = obj->p + obj->q + obj->M;
	if (obj->method == 0 || obj->method == 1) {
		printf("\n\n Exit Status \n");
		printf("Return Code : %d \n", obj->retval);
		printf("Exit Message : ");

		if (obj->retval == 0) {
			printf("Input Error");
		}
		else if (obj->retval == 1) {
			printf("Probable Success");
		}
		else if (obj->retval == 4) {
			printf("Optimization Routine didn't converge");
		}
		else if (obj->retval == 15) {
			printf("Optimization Routine Encountered Inf/Nan Values");
		}
	}
	printf("\n\n");
	printf(" ARIMA Order : ( %d, %d, %d) \n", obj->p, obj->d, obj->q);
	printf("\n");
	printf("%-20s%-20s%-20s \n\n", "Coefficients", "Value", "Standard Error");
	for (i = 0; i < obj->p; ++i) {
		printf("AR%-15d%-20g%-20g \n", i + 1, obj->phi[i],sqrt(obj->vcov[i+pq*i]));
	}
	for (i = 0; i < obj->q; ++i) {
		t = obj->p + i;
		printf("MA%-15d%-20g%-20g \n", i + 1, obj->theta[i], sqrt(obj->vcov[t + pq * t]));
	}
	printf("\n");
	t = obj->p + obj->q;
	if (obj->M == 1) {
		printf("%-17s%-20g%-20g \n", "MEAN", obj->mean, sqrt(obj->vcov[t + pq * t]));
	}
	else {
		printf("%-17s%-20g \n", "MEAN", obj->mean);
	}
	printf("\n");
	printf("%-17s%-20g \n", "SIGMA^2", obj->var);
	printf("\n");
	printf("ESTIMATION METHOD : ");
	if (obj->method == 0) {
		printf("MLE");
	}
	else if (obj->method == 1) {
		printf("CSS");
	}
	else if (obj->method == 2) {
		printf("BOX-JENKINS");
	}
	printf("\n\n");
	printf("OPTIMIZATION METHOD : ");
	if (obj->method == 2) {
		printf("Newton-Raphson");
	}
	else {
		if (obj->optmethod == 0) {
			printf("Nelder-Mead");
		}
		else if (obj->optmethod == 1) {
			printf("Newton Line Search");
		}
		else if (obj->optmethod == 2) {
			printf("Newton Trust Region - Hook Step");
		}
		else if (obj->optmethod == 3) {
			printf("Newton Trust Region - Double Dog-Leg");
		}
		else if (obj->optmethod == 4) {
			printf("Conjugate Gradient");
		}
		else if (obj->optmethod == 5) {
			printf("BFGS");
		}
		else if (obj->optmethod == 6) {
			printf("L-BFGS");
		}
		else if (obj->optmethod == 7) {
			printf("BFGS More-Thuente Line Search");
		}
	}
	printf("\n\n");
	printf("EQUATION FORM : x[t] ");
	for (i = 0; i < obj->p; ++i) {
		if (obj->phi[i] > 0) {
			printf("- %g*x[t - %d%s", fabs(obj->phi[i]), i + 1, "]");
		}
		else if (obj->phi[i] < 0) {
			printf("+ %g*x[t - %d%s", fabs(obj->phi[i]), i + 1, "] ");
		}
	}

	printf("=");
	if (obj->mean != 0.0) {
		printf(" %g +", obj->mean);
	}
	printf(" a[t] ");
	for (i = 0; i < obj->q; ++i) {
		if (obj->theta[i] > 0) {
			printf("- %g*a[t - %d%s", fabs(obj->theta[i]), i + 1, "]");
		}
		else if (obj->phi[i] < 0) {
			printf("+ %g*a[t - %d%s", fabs(obj->theta[i]), i + 1, "] ");
		}
	}
	printf("\n\n");
	if (obj->method == 0 || obj->method == 1) {
		printf("Log Likelihood : %g ", obj->loglik);
		printf("\n\n");
	}
	else {
		printf("Log Likelihood : Unavailable ");
		printf("\n\n");
	}
	if (obj->method == 0) {
		printf("AIC criterion : %g ", obj->aic);
		printf("\n\n");
	}
	else {
		printf("AIC Criterion : Unavailable ");
		printf("\n\n");
	}
}

void sarimax_summary(sarimax_object obj) {
	int i, pq,t,nd,ncxreg,mean;
	pq = obj->p + obj->q + obj->P + obj->Q + obj->M;
	mean = obj->M - obj->r;
	
	if (obj->method == 0 || obj->method == 1) {
		printf("\n\n Exit Status \n");
		printf("Return Code : %d \n", obj->retval);
		printf("Exit Message : ");

		if (obj->retval == 0) {
			printf("Input Error");
		}
		else if (obj->retval == 1) {
			printf("Probable Success");
		}
		else if (obj->retval == 4) {
			printf("Optimization Routine didn't converge");
		}
		else if (obj->retval == 7) {
			printf("Exogenous Variables are collinear");
		}
		else if (obj->retval == 10) {
			printf("Nonstationary AR part");
		}
		else if (obj->retval == 12) {
			printf("Nonstationary Seasonal AR part");
		}
		else if (obj->retval == 15) {
			printf("Optimization Routine Encountered Inf/Nan Values");
		}
	}
	printf("\n\n");
	printf("  ARIMA Seasonal Order : ( %d, %d, %d) * (%d, %d, %d) \n",obj->p,obj->d,obj->q, obj->P,obj->D,obj->Q );
	printf("\n");
	//mdisplay(obj->vcov,pq,pq);
	printf("%-20s%-20s%-20s \n\n", "Coefficients", "Value", "Standard Error");
	for (i = 0; i < obj->p; ++i) {
		printf("AR%-15d%-20g%-20g \n", i + 1, obj->phi[i], sqrt(obj->vcov[i + pq*i]));
	}
	for (i = 0; i < obj->q; ++i) {
		t = obj->p + i;
		printf("MA%-15d%-20g%-20g \n", i + 1, obj->theta[i], sqrt(obj->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->P; ++i) {
		t = obj->p + obj->q + i;
		printf("SAR%-14d%-20g%-20g \n", i + 1, obj->PHI[i], sqrt(obj->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->Q; ++i) {
		t = obj->p + obj->q + obj->P + i;
		printf("SMA%-14d%-20g%-20g \n", i + 1, obj->THETA[i], sqrt(obj->vcov[t + pq * t]));
	}
	printf("\n");
	t = obj->p + obj->q + obj->P + obj->Q;
	if (mean > 0) {
		printf("%-17s%-20g%-20g \n", "MEAN", obj->mean, sqrt(obj->vcov[t + pq * t]));
		for(i = 1; i <= obj->r; ++i) {
			t = obj->p + obj->q + obj->P + obj->Q + i;
			printf("%-17s%-20g%-20g \n", "EXOG", obj->exog[i-1], sqrt(obj->vcov[t + pq * t]));
		}
	}
	else {
		for(i = 0; i < obj->r; ++i) {
			t = obj->p + obj->q + obj->P + obj->Q + i;
			printf("%-17s%-20g%-20g \n", "EXOG", obj->exog[i], sqrt(obj->vcov[t + pq * t]));
		}
		printf("%-17s%-20g \n", "MEAN", obj->mean);
	}
	printf("\n");
	printf("%-17s%-20g \n", "SIGMA^2", obj->var);
	printf("\n");
	printf("ESTIMATION METHOD : ");
	if (obj->method == 0) {
		printf("CSS-MLE");
	}
	else if (obj->method == 1) {
		printf("MLE");
	}
	else if (obj->method == 2) {
		printf("CSS");
	}
	printf("\n\n");
	printf("OPTIMIZATION METHOD : ");
	
	if (obj->optmethod == 0) {
		printf("Nelder-Mead");
	}
	else if (obj->optmethod == 1) {
		printf("Newton Line Search");
	}
	else if (obj->optmethod == 2) {
		printf("Newton Trust Region - Hook Step");
	}
	else if (obj->optmethod == 3) {
		printf("Newton Trust Region - Double Dog-Leg");
	}
	else if (obj->optmethod == 4) {
		printf("Conjugate Gradient");
	}
	else if (obj->optmethod == 5) {
		printf("BFGS");
	}
	else if (obj->optmethod == 6) {
		printf("L-BFGS");
	}
	else if (obj->optmethod == 7) {
		printf("BFGS More-Thuente Line Search");
	}
	
	printf("\n\n");
	printf("EQUATION FORM : x[t] ");
	for (i = 0; i < obj->p; ++i) {
		if (obj->phi[i] > 0) {
			printf("- %g*x[t - %d%s", fabs(obj->phi[i]), i + 1, "]");
		}
		else if (obj->phi[i] < 0) {
			printf("+ %g*x[t - %d%s", fabs(obj->phi[i]), i + 1, "] ");
		}
	}

	printf("=");
	if (obj->mean != 0.0) {
		printf(" %g +", obj->mean);
	}
	printf(" a[t] ");
	for (i = 0; i < obj->q; ++i) {
		if (obj->theta[i] > 0) {
			printf("- %g*a[t - %d%s", fabs(obj->theta[i]), i + 1, "]");
		}
		else if (obj->phi[i] < 0) {
			printf("+ %g*a[t - %d%s", fabs(obj->theta[i]), i + 1, "] ");
		}
	}
	printf("\n\n");
	if (obj->method == 0 || obj->method == 1 || obj->method == 2) {
		printf("Log Likelihood : %g ", obj->loglik);
		printf("\n\n");
	}
	if (obj->method == 0 || obj->method == 1) {
		printf("AIC criterion : %g ", obj->aic);
		printf("\n\n");
	}
	else {
		printf("AIC Criterion : Unavailable ");
		printf("\n\n");
	}
}

void auto_arima_summary(auto_arima_object obj) {
	int i, pq,t,nd,ncxreg,mean,drift;
	pq = obj->p + obj->q + obj->P + obj->Q + obj->M;
	mean = obj->M - obj->r;

	if (obj->method == 0 || obj->method == 1) {
		printf("\n\n Exit Status \n");
		printf("Return Code : %d \n", obj->retval);
		printf("Exit Message : ");

		if (obj->retval == 0) {
			printf("Input Error");
		}
		else if (obj->retval == 1) {
			printf("Probable Success");
		}
		else if (obj->retval == 4) {
			printf("Optimization Routine didn't converge");
		}
		else if (obj->retval == 7) {
			printf("Exogenous Variables are collinear");
		}
		else if (obj->retval == 10) {
			printf("Nonstationary AR part");
		}
		else if (obj->retval == 12) {
			printf("Nonstationary Seasonal AR part");
		}
		else if (obj->retval == 15) {
			printf("Optimization Routine Encountered Inf/Nan Values");
		}
	}
	printf("\n\n");
	printf("  ARIMA Seasonal Order : ( %d, %d, %d) * (%d, %d, %d) \n",obj->p,obj->d,obj->q,
	 obj->P,obj->D,obj->Q );
	printf("\n");

	printf("%-20s%-20s%-20s \n\n", "Coefficients", "Value", "Standard Error");
	for (i = 0; i < obj->p; ++i) {
		printf("AR%-15d%-20g%-20g \n", i + 1, obj->phi[i], sqrt(obj->vcov[i + pq*i]));
	}
	for (i = 0; i < obj->q; ++i) {
		t = obj->p + i;
		printf("MA%-15d%-20g%-20g \n", i + 1, obj->theta[i], sqrt(obj->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->P; ++i) {
		t = obj->p + obj->q + i;
		printf("SAR%-14d%-20g%-20g \n", i + 1, obj->PHI[i], sqrt(obj->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->Q; ++i) {
		t = obj->p + obj->q + obj->P + i;
		printf("SMA%-14d%-20g%-20g \n", i + 1, obj->THETA[i], sqrt(obj->vcov[t + pq * t]));
	}
	printf("\n");
	t = obj->p + obj->q + obj->P + obj->Q;
	if (mean > 0) {
		printf("%-17s%-20g%-20g \n", "MEAN", obj->mean, sqrt(obj->vcov[t + pq * t]));
		t++;
		
	}
	else {
		printf("%-17s%-20g \n", "MEAN", obj->mean);
	}
	ncxreg = 0;
	if (obj->idrift == 1) {
		printf("%-17s%-20g%-20g \n", "TREND", obj->exog[0], sqrt(obj->vcov[t + pq * t]));
		t++;
		ncxreg++;
	} else {
		printf("%-17s%-20g \n", "TREND", 0.0);
	}
	for(i = ncxreg; i < obj->r; ++i) {
		printf("%-17s%-20g%-20g \n", "EXOG", obj->exog[i], sqrt(obj->vcov[t + pq * t]));
		t++;
	}
	printf("\n");
	printf("%-17s%-20g \n", "SIGMA^2", obj->sigma2);
	printf("\n");
	printf("ESTIMATION METHOD : ");
	if (obj->method == 0) {
		printf("CSS-MLE");
	}
	else if (obj->method == 1) {
		printf("MLE");
	}
	else if (obj->method == 2) {
		printf("CSS");
	}
	printf("\n\n");
	printf("OPTIMIZATION METHOD : ");
	
	if (obj->optmethod == 0) {
		printf("Nelder-Mead");
	}
	else if (obj->optmethod == 1) {
		printf("Newton Line Search");
	}
	else if (obj->optmethod == 2) {
		printf("Newton Trust Region - Hook Step");
	}
	else if (obj->optmethod == 3) {
		printf("Newton Trust Region - Double Dog-Leg");
	}
	else if (obj->optmethod == 4) {
		printf("Conjugate Gradient");
	}
	else if (obj->optmethod == 5) {
		printf("BFGS");
	}
	else if (obj->optmethod == 6) {
		printf("L-BFGS");
	}
	else if (obj->optmethod == 7) {
		printf("BFGS More-Thuente Line Search");
	}
	
	printf("\n\n");
	
	printf("AIC criterion : %g ", obj->aic);
	printf("\n\n");
	printf("BIC criterion : %g ", obj->bic);
	printf("\n\n");
	printf("AICC criterion : %g ", obj->aicc);
	printf("\n\n");
	if (obj->method == 0 || obj->method == 1 || obj->method == 2) {
		printf("Log Likelihood : %g ", obj->loglik);
		printf("\n\n");
	}
}

void sarimax_wrapper_summary(sarimax_wrapper_object obj) {
	int i, pq,t,nd,ncxreg,mean;
	pq = obj->sarimax->p + obj->sarimax->q + obj->sarimax->P + obj->sarimax->Q + obj->sarimax->M;
	mean = obj->sarimax->M - obj->sarimax->r;
	
	if (obj->sarimax->method == 0 || obj->sarimax->method == 1) {
		printf("\n\n Exit Status \n");
		printf("Return Code : %d \n", obj->sarimax->retval);
		printf("Exit Message : ");

		if (obj->sarimax->retval == 0) {
			printf("Input Error");
		}
		else if (obj->sarimax->retval == 1) {
			printf("Probable Success");
		}
		else if (obj->sarimax->retval == 4) {
			printf("Optimization Routine didn't converge");
		}
		else if (obj->sarimax->retval == 7) {
			printf("Exogenous Variables are collinear");
		}
		else if (obj->sarimax->retval == 10) {
			printf("Nonstationary AR part");
		}
		else if (obj->sarimax->retval == 12) {
			printf("Nonstationary Seasonal AR part");
		}
		else if (obj->sarimax->retval == 15) {
			printf("Optimization Routine Encountered Inf/Nan Values");
		}
	}
	printf("\n\n");
	printf("  ARIMA Seasonal Order : ( %d, %d, %d) * (%d, %d, %d) \n",obj->sarimax->p,obj->sarimax->d,obj->sarimax->q,
	 obj->sarimax->P,obj->sarimax->D,obj->sarimax->Q );
	printf("\n");
	//mdisplay(obj->vcov,pq,pq);
	printf("%-20s%-20s%-20s \n\n", "Coefficients", "Value", "Standard Error");
	for (i = 0; i < obj->sarimax->p; ++i) {
		printf("AR%-15d%-20g%-20g \n", i + 1, obj->sarimax->phi[i], sqrt(obj->sarimax->vcov[i + pq*i]));
	}
	for (i = 0; i < obj->sarimax->q; ++i) {
		t = obj->sarimax->p + i;
		printf("MA%-15d%-20g%-20g \n", i + 1, obj->sarimax->theta[i], sqrt(obj->sarimax->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->sarimax->P; ++i) {
		t = obj->sarimax->p + obj->sarimax->q + i;
		printf("SAR%-14d%-20g%-20g \n", i + 1, obj->sarimax->PHI[i], sqrt(obj->sarimax->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->sarimax->Q; ++i) {
		t = obj->sarimax->p + obj->sarimax->q + obj->sarimax->P + i;
		printf("SMA%-14d%-20g%-20g \n", i + 1, obj->sarimax->THETA[i], sqrt(obj->sarimax->vcov[t + pq * t]));
	}
	printf("\n");
	t = obj->sarimax->p + obj->sarimax->q + obj->sarimax->P + obj->sarimax->Q;
	if (mean > 0) {
		printf("%-17s%-20g%-20g \n", "MEAN", obj->sarimax->mean, sqrt(obj->sarimax->vcov[t + pq * t]));
		t++;
		
	}
	else {
		printf("%-17s%-20g \n", "MEAN", obj->sarimax->mean);
	}
	ncxreg = 0;
	if (obj->idrift == 1) {
		printf("%-17s%-20g%-20g \n", "TREND", obj->sarimax->exog[0], sqrt(obj->sarimax->vcov[t + pq * t]));
		t++;
		ncxreg++;
	} else {
		printf("%-17s%-20g \n", "TREND", 0.0);
	}
	for(i = ncxreg; i < obj->sarimax->r; ++i) {
		printf("%-17s%-20g%-20g \n", "EXOG", obj->sarimax->exog[i], sqrt(obj->sarimax->vcov[t + pq * t]));
		t++;
	}
	printf("\n");
	printf("%-17s%-20g \n", "SIGMA^2", obj->sigma2);
	printf("\n");
	printf("ESTIMATION METHOD : ");
	if (obj->sarimax->method == 0) {
		printf("CSS-MLE");
	}
	else if (obj->sarimax->method == 1) {
		printf("MLE");
	}
	else if (obj->sarimax->method == 2) {
		printf("CSS");
	}
	printf("\n\n");
	printf("OPTIMIZATION METHOD : ");
	
	if (obj->sarimax->optmethod == 0) {
		printf("Nelder-Mead");
	}
	else if (obj->sarimax->optmethod == 1) {
		printf("Newton Line Search");
	}
	else if (obj->sarimax->optmethod == 2) {
		printf("Newton Trust Region - Hook Step");
	}
	else if (obj->sarimax->optmethod == 3) {
		printf("Newton Trust Region - Double Dog-Leg");
	}
	else if (obj->sarimax->optmethod == 4) {
		printf("Conjugate Gradient");
	}
	else if (obj->sarimax->optmethod == 5) {
		printf("BFGS");
	}
	else if (obj->sarimax->optmethod == 6) {
		printf("L-BFGS");
	}
	else if (obj->sarimax->optmethod == 7) {
		printf("BFGS More-Thuente Line Search");
	}
	
	printf("\n\n");
	
	printf("AIC criterion : %g ", obj->aic);
	printf("\n\n");
	printf("BIC criterion : %g ", obj->bic);
	printf("\n\n");
	printf("AICC criterion : %g ", obj->aicc);
	printf("\n\n");
	if (obj->sarimax->method == 0 || obj->sarimax->method == 1 || obj->sarimax->method == 2) {
		printf("Log Likelihood : %g ", obj->sarimax->loglik);
		printf("\n\n");
	}
}

void myarima_summary(myarima_object obj) {
	int i, pq,t,nd,ncxreg,mean;
	pq = obj->sarimax->p + obj->sarimax->q + obj->sarimax->P + obj->sarimax->Q + obj->sarimax->M;
	mean = obj->sarimax->M - obj->sarimax->r;
	
	if (obj->sarimax->method == 0 || obj->sarimax->method == 1) {
		printf("\n\n Exit Status \n");
		printf("Return Code : %d \n", obj->sarimax->retval);
		printf("Exit Message : ");

		if (obj->sarimax->retval == 0) {
			printf("Input Error");
		}
		else if (obj->sarimax->retval == 1) {
			printf("Probable Success");
		}
		else if (obj->sarimax->retval == 4) {
			printf("Optimization Routine didn't converge");
		}
		else if (obj->sarimax->retval == 7) {
			printf("Exogenous Variables are collinear");
		}
		else if (obj->sarimax->retval == 10) {
			printf("Nonstationary AR part");
		}
		else if (obj->sarimax->retval == 12) {
			printf("Nonstationary Seasonal AR part");
		}
		else if (obj->sarimax->retval == 15) {
			printf("Optimization Routine Encountered Inf/Nan Values");
		}
	}
	printf("\n\n");
	printf("  ARIMA Seasonal Order : ( %d, %d, %d) * (%d, %d, %d) \n",obj->sarimax->p,obj->sarimax->d,obj->sarimax->q,
	 obj->sarimax->P,obj->sarimax->D,obj->sarimax->Q );
	printf("\n");
	//mdisplay(obj->vcov,pq,pq);
	printf("%-20s%-20s%-20s \n\n", "Coefficients", "Value", "Standard Error");
	for (i = 0; i < obj->sarimax->p; ++i) {
		printf("AR%-15d%-20g%-20g \n", i + 1, obj->sarimax->phi[i], sqrt(obj->sarimax->vcov[i + pq*i]));
	}
	for (i = 0; i < obj->sarimax->q; ++i) {
		t = obj->sarimax->p + i;
		printf("MA%-15d%-20g%-20g \n", i + 1, obj->sarimax->theta[i], sqrt(obj->sarimax->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->sarimax->P; ++i) {
		t = obj->sarimax->p + obj->sarimax->q + i;
		printf("SAR%-14d%-20g%-20g \n", i + 1, obj->sarimax->PHI[i], sqrt(obj->sarimax->vcov[t + pq * t]));
	}
	for (i = 0; i < obj->sarimax->Q; ++i) {
		t = obj->sarimax->p + obj->sarimax->q + obj->sarimax->P + i;
		printf("SMA%-14d%-20g%-20g \n", i + 1, obj->sarimax->THETA[i], sqrt(obj->sarimax->vcov[t + pq * t]));
	}
	printf("\n");
	t = obj->sarimax->p + obj->sarimax->q + obj->sarimax->P + obj->sarimax->Q;
	if (mean > 0) {
		printf("%-17s%-20g%-20g \n", "MEAN", obj->sarimax->mean, sqrt(obj->sarimax->vcov[t + pq * t]));
		t++;
		
	}
	else {
		printf("%-17s%-20g \n", "MEAN", obj->sarimax->mean);
	}
	ncxreg = 0;
	if (obj->idrift == 1) {
		printf("%-17s%-20g%-20g \n", "TREND", obj->sarimax->exog[0], sqrt(obj->sarimax->vcov[t + pq * t]));
		t++;
		ncxreg++;
	} else {
		printf("%-17s%-20g \n", "TREND", 0.0);
	}
	for(i = ncxreg; i < obj->sarimax->r; ++i) {
		printf("%-17s%-20g%-20g \n", "EXOG", obj->sarimax->exog[i], sqrt(obj->sarimax->vcov[t + pq * t]));
		t++;
	}
	printf("\n");
	printf("%-17s%-20g \n", "SIGMA^2", obj->sigma2);
	printf("\n");
	printf("ESTIMATION METHOD : ");
	if (obj->sarimax->method == 0) {
		printf("CSS-MLE");
	}
	else if (obj->sarimax->method == 1) {
		printf("MLE");
	}
	else if (obj->sarimax->method == 2) {
		printf("CSS");
	}
	printf("\n\n");
	printf("OPTIMIZATION METHOD : ");
	
	if (obj->sarimax->optmethod == 0) {
		printf("Nelder-Mead");
	}
	else if (obj->sarimax->optmethod == 1) {
		printf("Newton Line Search");
	}
	else if (obj->sarimax->optmethod == 2) {
		printf("Newton Trust Region - Hook Step");
	}
	else if (obj->sarimax->optmethod == 3) {
		printf("Newton Trust Region - Double Dog-Leg");
	}
	else if (obj->sarimax->optmethod == 4) {
		printf("Conjugate Gradient");
	}
	else if (obj->sarimax->optmethod == 5) {
		printf("BFGS");
	}
	else if (obj->sarimax->optmethod == 6) {
		printf("L-BFGS");
	}
	else if (obj->sarimax->optmethod == 7) {
		printf("BFGS More-Thuente Line Search");
	}
	
	printf("\n\n");
	
	printf("AIC criterion : %g ", obj->aic);
	printf("\n\n");
	printf("BIC criterion : %g ", obj->bic);
	printf("\n\n");
	printf("AICC criterion : %g ", obj->aicc);
	printf("\n\n");
	if (obj->sarimax->method == 0 || obj->sarimax->method == 1 || obj->sarimax->method == 2) {
		printf("Log Likelihood : %g ", obj->sarimax->loglik);
		printf("\n\n");
	}
}

void aa_ret_summary(aa_ret_object obj) {
	if (obj->otype == 2) {
		sarimax_wrapper_summary(obj->Arima);
	} else if (obj->otype == 1) {
		myarima_summary(obj->myarima);
	} else {
		printf("aa_ret Error \n");
	}
}

int ar_estimate(double *x, int N, int method) {
	int p, ordermax, logn, i;
	double wmean, aic, loglik, var, lvar, aic0,cssml;
	double *inp, *resid, *hess, *phi;

	inp = (double*)malloc(sizeof(double)* N);
	resid = (double*)malloc(sizeof(double)* N);

	logn = (int) (10.0 * log((double)N));
	ordermax = imin(N - 1, logn);

	if (method == 2) {
		ordermax = imin(ordermax, 12); // MLE
	}
	wmean = mean(x, N);
	aic0 = 0.0;
	cssml = 0;

	hess = (double*)malloc(sizeof(double)* (ordermax + 1) * (ordermax + 1));
	phi = (double*)malloc(sizeof(double)* ordermax);

	for (i = 0; i < N; ++i) {
		inp[i] = x[i] - wmean;
	}

	for (i = 0; i < ordermax; ++i) {
		if (method == 0) {
			ywalg2(inp, N, i + 1, phi, &var);
		}
		else if (method == 1) {
			burgalg(inp, N-1, i + 1, phi, &var);
		}
		else if (method == 2) {
			as154(inp, N, 7, i + 1, 0, 0, phi, NULL, &wmean, &var, resid, &loglik, hess,cssml);
		}
		else {
			printf("\n The code only accepts numerical values 0,1 and 2 \n");
			printf("\n Method 0 : Yule-Walker \n");
			printf("\n Method 1 : Burg \n");
			printf("\n Method 2 : MLE \n");
		}
		lvar = log(var);
		aic = lvar + 2 * (double)(i + 1) / N;
		if (i == 0) {
			aic = aic0;
			p = 1;
		}
		else {
			if (aic < aic0) {
				aic0 = aic;
				p = i + 1;
			}
		}
	}

	printf("\n\n");
	printf("AIC Estimate : p = %d \n", p);

	free(inp);
	free(resid);
	free(hess);
	free(phi);
	return p;
}

void ar_summary(ar_object obj) {
	int i;

	printf("\n\n");
	printf(" AR Order : ( %d) \n", obj->p);
	printf("\n");
	printf("%-20s%-20s \n\n", "Coefficients", "Value");
	for (i = 0; i < obj->p; ++i) {
		printf("AR%-15d%-20g \n", i + 1, obj->phi[i]);
	}
	printf("\n");
	
	printf("%-17s%-20g \n", "MEAN", obj->mean);
	
	printf("\n");
	printf("%-17s%-20g \n", "SIGMA^2", obj->var);
	printf("\n");
	printf("ESTIMATION METHOD : ");
	if (obj->method == 0) {
		printf("Yule Walker");
	}
	else if (obj->method == 1) {
		printf("Burg");
	}
	else if (obj->method == 2) {
		printf("MLE");
	}
	printf("\n\n");
}

void model_estimate(double *x, int N, int d, int pmax, int h) {
	int m, i, t, pq, j;
	double wmean, sos, var, lvar, aic, sc, hq, aic0, sc0, hq0;
	double *inp, *phim, *a;
	int paic, qaic, psc, qsc, phq, qhq;
	/*
	x - Input Time Series of length N
	d - Number of time the series has to be differenced (d = 0 when there are no trends)
	pmax - Maximum AR and MA order pmax = max(p,q)
	h - Order of AR model to be fitted to obtain residuals
	0 <= pmax <= h
	Typically pmax = 3,4 and h = 8. Increase value of h if you want to fit high order ARMA process
	*/

	inp = (double*)malloc(sizeof(double)* (N - d));

	if (d > 0) {
		N = diff(x, N, d, inp); // No need to demean x
	}

	m = h;

	phim = (double*)malloc(sizeof(double)* m);
	a = (double*)malloc(sizeof(double)* N);

	wmean = mean(x, N);

	for (i = 0; i < N; ++i) {
		inp[i] = x[i] - wmean;
	}

	// Estimate AR(m) coefficients

	ywalg(inp, N, m, phim);

	for (i = 0; i < N; ++i) {
		a[i] = 0.0;
	}

	for (t = m; t < N; ++t) {
		a[t] = inp[t];
		for (i = 0; i < m; ++i) {
			a[t] -= phim[i] * inp[t - i - 1];
		}
	}
	aic0 = sc0 = hq0 = 0.0;
	paic = qaic = psc = qsc = phq = qhq = 0;

	for (i = 0; i <= pmax; ++i) {
		for (j = 0; j <= pmax; ++j) {
			pq = i + j;
			if (pq >= 1) {
				hrstep2(inp, N, m, i, j, pq, a, &sos, &var);
				lvar = log(var);
				aic = lvar + 2 * (double)pq / N;
				sc = lvar + log((double)N) * (double)pq / N;
				hq = lvar + 2 * log(log((double)N)) * (double)pq / N;
				if (i == 0 && j == 1) {
					aic0 = aic;
					sc0 = sc;
					hq0 = hq;
					paic = psc = phq = 0;
					qaic = qsc = qhq = 1;
				}
				else {
					if (aic < aic0) {
						aic0 = aic;
						paic = i;
						qaic = j;
					}
					if (sc < sc0) {
						sc0 = sc;
						psc = i;
						qsc = j;
					}
					if (hq < hq0) {
						hq0 = hq;
						phq = i;
						qhq = j;
					}
				}
				printf("\n p = %d q = %d AIC = %g SC = %g HQ = %g \n", i, j, aic, sc, hq);
			}
		}
	}

	printf("\n\n");
	printf("AIC Estimate : p = %d q = %d \n", paic, qaic);
	printf("SC Estimate  : p = %d q = %d \n", psc, qsc);
	printf("HQ Estimate  : p = %d q = %d \n", phq, qhq);

	free(inp);
	free(phim);
	free(a);
}

void pacf_opt(double* vec, int N, int method, double* par, int M) {
	// Method 0 : Yule-Walker
	// Method 1 : Burg
	// Method 2 : Box-Jenkins Conditional MLE

	if (method == 0) {
		pacf_yw(vec, N, par, M);
	}
	else if (method == 1) {
		pacf_burg(vec, N, par, M);
	}
	else if (method == 2) {
		pacf_mle(vec, N, par, M);
	}
	else {
		printf("\n The code only accepts numerical values 0,1 and 2 \n");
		printf("\n Method 0 : Yule-Walker \n");
		printf("\n Method 1 : Burg \n");
		printf("\n Method 2 : Box-Jenkins Conditional MLE \n");
	}
}

void pacf(double* vec, int N, double* par, int M) {
	pacf_yw(vec, N, par, M);
	/*
	auto_fft_object obj;
	double *acorr;

	acorr = (double*)malloc(sizeof(double)* M);

	if (N < 50) {
	autocorr(vec, N, acorr, M);
	}
	else {
	obj = auto_fft_init(N);
	autocorr_fft(obj, vec, acorr, M);
	free_auto(obj);
	}

	free(acorr);
	*/
}

void acvf(double* vec, int N, double* par, int M) {
	// Auto Covariance Function
	autocovar(vec, N, par, M);
}

void acvf_opt(double* vec, int N, int method, double* par, int M) {
	//Method 0 : Regular Method. Slow for large input length.
	//Method 1 : FFT based Method. Use it if data length is large
	auto_fft_object obj;

	if (method == 0) {
		autocovar(vec, N, par, M);
	}
	else if (method == 1) {
		obj = auto_fft_init(N);
		autocovar_fft(obj, vec, par, M);
		free_auto(obj);
	}
	else {
		printf("\n The code only accepts numerical values 0 and 1 \n");
		printf("\n Method 0 : General Method \n");
		printf("\n Method 1 : FFT \n");

	}
}

void acvf2acf(double *acorr, int M) {
	// Converts Autocovariance to autocorrelation function
	int i;
	double var;

	var = acorr[0];
	acorr[0] = 1.0;

	for (i = 1; i < M; i++) {
		acorr[i] = acorr[i] / var;
	}
}

void arima_free(arima_object object) {
	free(object);
}

void sarimax_free(sarimax_object object) {
	free(object);
}

void sarimax_wrapper_free(sarimax_wrapper_object object) {
	sarimax_free(object->sarimax);
	free(object);
}

void myarima_free(myarima_object object) {
	sarimax_free(object->sarimax);
	free(object);
}

void aa_ret_free(aa_ret_object object) {
	if (object->otype == 1) {
		free(object->myarima);
	} else {
		free(object->Arima);
	}
	free(object);
}

void sarima_free(sarima_object object) {
	free(object);
}

void ar_free(ar_object object) {
	free(object);
}

// Yule-Walker, Burg and Hannan Rissanen Algorithms for Initial Parameter Estimation

void yw(double *x, int N, int p, double *phi, double *var) {
	ywalg2(x, N, p, phi, var);
}

void burg(double *x, int N, int p, double *phi, double *var) {
	burgalg(x, N-1, p, phi, var);
}

void hr(double *x, int N, int p, int q, double *phi, double *theta, double *var) {
	hralg(x, N, p, q, phi, theta, var);
}