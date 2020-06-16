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

sarimax_object sarimax_init(int p, int d, int q,int P, int D, int Q,int s, int r, int N) {
	sarimax_object obj = NULL;
	int i,ncxreg;
	if (P + D + Q == 0) {
		s = 0;
	}

	if (d + D > 0) {
		ncxreg = r;
	}
	else {
		ncxreg = 1 + r;
	}
	if (p < 0 || d < 0 || q < 0 || N <= 0 || P < 0 || D < 0 || Q < 0) {
		printf("\n Input Values cannot be Negative. Program Exiting. \n");
		exit(-1);
	}
	//retval = 0 Input Error
	// retval = 1 Probable Success
	// retval = 4 Optimization Routine didn't converge
	// Retval = 15 Optimization Routine Encountered Inf/Nan Values
	
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
	int p,q,d,N,M,P,D,Q,s,r,nd,ncxreg;
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
		obj->retval = as154x(inp, obj->N, xreg, obj->optmethod, obj->p, obj->d, obj->q, obj->s, obj->P, obj->D, obj->Q, obj->params, obj->params + p, obj->params + p + q, 
			obj->params + p + q + P,obj->params + p + q + P + Q, obj->r, &obj->mean, &obj->var,obj->params + p + q + P + Q + M,
			 &obj->loglik, obj->params + p + q + P + Q + M + N - d - s*D);
		//mdisplay(obj->vcov,p+q+P+Q+M,p+q+P+Q+M);
		obj->loglik = -0.5 * (obj->Nused * (2 * obj->loglik + 1.0 + log(2 * 3.14159)));
		obj->aic = -2.0 * obj->loglik + 2.0 * (obj->p + obj->q + obj->P + obj->Q + obj->M) + 2.0;
	}
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
	//mdisplay(delta, 1, d + D*s + 1);
	
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