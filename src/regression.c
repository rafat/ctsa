// SPDX-License-Identifier: BSD-3-Clause
/*
 * regression.c
 *
 *  Created on: Jun 5, 2013
 *      Author: Rafat Hussain
 */

#include "regression.h"


reg_object reg_init(int N, int p) {
	reg_object obj = NULL;
	int i;
	/*
	* reg_init initializes the regression object.
	* N - number of observation samples. All variables (indepnedent and dependent)
	* are of identical length N.
	*
	* p - Total Number of variables (independent and dependent).
	* The base case is p = 2 corresonding to y = b0 + b1 *x + u
	*/
	if (p < 0) {
		printf("The base case is p = 1 corresonding to either y = b0 + u or y = b1 *x + u \n");
		printf("p = Number of Independent Variables (p-1) + Dependent Variable (1) \n");
		exit(1);
	}
	obj = (reg_object)malloc(sizeof(struct reg_set) + sizeof(bparam)* (p + 1));
	obj->N = N;
	obj->p = p;
	obj->alpha = 0.05;
	obj->sigma = 0.0;
	obj->sigma_lower = 0.0;
	obj->sigma_upper = 0.0;
	obj->r2 = obj->R2[0] = 0.0;
	obj->r2adj = obj->R2[1] = 0.0;
	obj->df = 0;
	obj->intercept = (p == 0) ? 0 : 1;
	strcpy(obj->lls, "qr");

	obj->TSS = 0.0;
	obj->ESS = 0.0;
	obj->RSS = 0.0;
	obj->df_ESS = 0;
	obj->df_RSS = 0;
	obj->FStat = 0.0;
	obj->PVal = 0.0;
	obj->loglik = 0.0;
	obj->aic = 0.0;
	obj->bic = 0.0;
	obj->aicc = 0.0;

	for (i = 0; i < p + 1; ++i) {
		(obj->beta + i)->value = 0.0;
		(obj->beta + i)->lower = 0.0;
		(obj->beta + i)->upper = 0.0;
		(obj->beta + i)->stdErr = 0.0;
	}


	return obj;
}

void linreg_clrm(double *x, double *y, int N, double* b,
	double *var, double *res, double alpha, double *anv,
	double* ci_lower, double* ci_upper) {
	int i, df;
	double *xi;
	double *yi;
	double xm, ym, num, den, sum, den2, den3;
	double seb1, seb2, a2, ta2, intrvl, chia2, chi1a2;

	/*
	* Classic Linear Regression Model
	*
	* x,y and res are all of same length - N
	*
	* alpha is used to determine confidence interval limits
	* for a given 100*(1-alpha) % confidence interval
	*
	*
	* Alpha takes values between  0 and 1.
	*
	* For 95% confidence interval, the value of alpha is 0.05
	* For 90% confidence interval, the value of alpha is 0.10 etc.
	*/

	/*
	* Parameters
	* b[0] - beta1
	* b[1] - beta2
	*/

	/*
	* var[0] - variance of residuals
	* var[1] variance beta1
	* var[2] variance beta2
	* var[3] covariance beta1,beta2
	* var[4] r^2 Goodness of Fit
	*/

	/*
	* ANOVA
	*
	* anv[0] - TSS Total Sum Of Squares
	* anv[1] - ESS Explained Sum Of Squares
	* anv[2] - RSS Residual Sum Of Squares
	* anv[3] - degrees of freedom of ESS
	* anv[4] - degrees of freedom of RSS
	* anv[5] - F Statistics = (anv[1] / anv[3]) / (anv[2] / anv[4])
	* anv[6] - P value associated with anv[5] used to reject/accept
	* zero hypothesis
	*/

	xi = (double*)malloc(sizeof(double) * N);
	yi = (double*)malloc(sizeof(double) * N);

	xm = mean(x, N);
	ym = mean(y, N);
	num = 0.0;
	den = 0.0;
	den2 = 0.0;
	den3 = 0.0;
	sum = 0.0;

	for (i = 0; i < N; ++i) {
		xi[i] = x[i] - xm;
		yi[i] = y[i] - ym;
		num += xi[i] * yi[i];
		den += xi[i] * xi[i];
		den2 += x[i] * x[i];
		den3 += yi[i] * yi[i];
	}

	b[1] = num / den;
	b[0] = ym - b[1] * xm;
	num = 0.0;

	for (i = 0; i < N; ++i) {
		res[i] = y[i] - b[0] - b[1] * x[i];
		sum += res[i] * res[i];
	}


	var[0] = sum / (N - 2);
	var[1] = den2 * var[0] / (N * den);
	var[2] = var[0] / den;
	var[3] = -xm * var[2];
	var[4] = b[1] * b[1] * den / den3;
	//Confidence Interval Estimations for a given alpha value

	// Find standard errors seb1 and seb2 from variances

	seb1 = sqrt(var[1]);
	seb2 = sqrt(var[2]);

	a2 = alpha / 2.;
	df = N - 2; // Number of parameters = 2

	ta2 = tinv(a2, df);
	intrvl = ta2 * seb1;

	ci_lower[0] = b[0] - intrvl;
	ci_upper[0] = b[0] + intrvl;

	intrvl = ta2 * seb2;

	ci_lower[1] = b[1] - intrvl;
	ci_upper[1] = b[1] + intrvl;

	// confidence interval for sigma^2

	chia2 = chiinv(a2, df);
	chi1a2 = chiinv(1. - a2, df);

	ci_lower[2] = (double)df * var[0] / chi1a2;
	ci_upper[2] = (double)df * var[0] / chia2;

	// ANOVA

	df = N - 2;
	anv[0] = den3;
	anv[2] = var[0] * (double)df;
	anv[1] = den3 - anv[2];
	anv[3] = 1;
	anv[4] = (double)df;
	anv[5] = (anv[1] / anv[3]) / (anv[2] / anv[4]);
	printf("%lf \n", var[0]);
	anv[6] = 1.0 - fcdf(anv[5], 1, df);

	free(xi);
	free(yi);
}

void zerohyp_clrm(int N, double *b, double *val, double *tval, double *pval) {
	int df;
	/*
	* Zero Hypothesis Test for CLRM
	*/

	df = N - 2;

	tval[0] = fabs(b[0]) / sqrt(val[1]);
	tval[1] = fabs(b[1]) / sqrt(val[2]);

	pval[0] = 1.0 - tcdf(fabs(tval[0]), df);
	pval[1] = 1.0 - tcdf(fabs(tval[1]), df);
}

void linreg_multi(int p, double *xi, double *y, int N, double* b, double *sigma2,
	double *xxti, double *R2, double *res, double alpha, double *anv,
	double* ci_lower, double* ci_upper, int *rank, char *llsmethod, int intercept) {
	/*
	p corresponds to number of coefficients (including the intercept)
	intercept - 1 : include intercept term. 0 : No intercept term
	*/
	double *x, *xt, *seb, *xb;
	double *xxt, *yxt, *temp1, *temp2, *bt, *xxt2;
	int i, df;
	int *ipiv;
	double tmp1, sum, ym, ym2, a2, ta2, intrvl, chia2, chi1a2, dfi, mss, rss;

	/*
	* Matrix Based Multiple Regression
	* y has the length N
	* xi has the length (p-1) X N. It contains p-1 independent variables
	* of length N each
	*/

	x = (double*)malloc(sizeof(double) * p * N);
	xt = (double*)malloc(sizeof(double) * p * N);
	xb = (double*)malloc(sizeof(double) * N);
	xxt = (double*)malloc(sizeof(double) * p * p);
	xxt2 = (double*)malloc(sizeof(double) * p * p);
	yxt = (double*)malloc(sizeof(double) * p);
	ipiv = (int*)malloc(sizeof(int) * p);
	//yt = (double*) malloc (sizeof(double) * N);
	temp1 = (double*)calloc(1, sizeof(double));
	temp2 = (double*)malloc(sizeof(double) * p);
	seb = (double*)malloc(sizeof(double) * p);
	bt = (double*)calloc(1, sizeof(double));

	// Initalize xt matrix - Transpose Matrix is initilaized first to
	// deal with column major nature of the X matrix

	if (intercept == 1) {
		for (i = 0; i < N; ++i) {
			xt[i] = 1.;
		}

		for (i = N; i < p*N; ++i) {
			xt[i] = xi[i - N];
		}
	}
	else if (intercept == 0) {
		for (i = 0; i < p*N; ++i) {
			xt[i] = xi[i];
		}
	}


	// Calculate x * x'

	mtranspose(xt, p, N, x);

	//cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ra,cb,ca,1.0,x,ca,xt,cb,0.0,xxt,cb);
	mmult(xt, x, xxt, p, N, p);
	memcpy(xxt2, xxt, sizeof(double)*p*p);

	// Calculate  x' * y

	//cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ra,cb,ca,1.0,y,ca,xt,cb,0.0,yxt,cb);
	mmult(xt, y, yxt, p, N, 1);

	// Find LU decomposition and then solve Ax = B
	// A = L * U where A = xxt
	// x is b vector and
	// B is yxt

	//ludecomp(xxt,p,ipiv);
	//linsolve(xxt,p,yxt,ipiv,b);
	if (p > 0) {
		if (!strcmp(llsmethod, "qr")) {
			*rank = lls_qr(xxt, yxt, p, p, b);
		}
		else if (!strcmp(llsmethod, "normal")) {
			lls_normal(xxt, yxt, p, p, b);
			printf("Warning - This method will not calculate the rank of the regression. Use method qr instead. \n");
		}
		else if (!strcmp(llsmethod, "svd")) {
			lls_svd2(xxt, yxt, p, p, b);
			printf("Warning - This method will not calculate the rank of the regression. Use method qr instead. \n");
		}
		else {
			printf("This function only accepts one of three least square methods - \n");
			printf(" qr, normal and svd. \n");
			exit(-1);
		}
	}

	//Find Variance of the residuals

	//mtranspose(y,N,1,yt);
	df = N - p;
	mmult(y, y, bt, 1, N, 1);
	mmult(xt, y, temp2, p, N, 1);
	mmult(b, temp2, temp1, 1, p, 1);

	*sigma2 = *bt - *temp1;
	*sigma2 /= (double)df;
	tmp1 = *sigma2;

	// Variance-Covariance Matrix of beta parameters
	ludecomp(xxt2, p, ipiv);
	minverse(xxt2, p, ipiv, xxti);

	for (i = 0; i < p*p; ++i) {
		xxti[i] *= tmp1;
	}


	// r^2 Goodness of fit
	sum = 0.;
	for (i = 0; i < N; ++i) {
		sum += y[i];
	}
	ym = sum / N;
	ym2 = ym * ym;

	R2[0] = (*temp1 - N*ym2) / (*bt - N*ym2);

	// Adjusted r^2

	dfi = intercept == 1 ? 1.0 : 0.0;

	R2[1] = 1.0 - (1.0 - R2[0]) * ((double)N - dfi) / ((double)df);

	a2 = alpha / 2.;

	//Confidence Interval for beta parameters

	for (i = 0; i < p; ++i) {
		seb[i] = sqrt(xxti[i*(p + 1)]);

		ta2 = tinv(a2, df);
		intrvl = ta2 * seb[i];

		ci_lower[i] = b[i] - intrvl;
		ci_upper[i] = b[i] + intrvl;

	}

	// confidence interval for sigma^2

	chia2 = chiinv(a2, df);
	chi1a2 = chiinv(1. - a2, df);

	ci_lower[p] = (double)df * tmp1 / chi1a2;
	ci_upper[p] = (double)df * tmp1 / chia2;

	//Residual
	/* u = y - Xb
	*/
	mmult(x, b, xb, N, p, 1);

	for (i = 0; i < N; ++i) {
		res[i] = y[i] - xb[i];
	}


	mss = 0.0;
	rss = 0.0;

	for (i = 0; i < N; ++i) {
		mss += (xb[i] * xb[i]);
		rss += (res[i] * res[i]);
	}

	// ANOVA

	/*
	* ANOVA
	*
	* anv[0] - TSS Total Sum Of Squares
	* anv[1] - ESS Explained Sum Of Squares
	* anv[2] - RSS Residual Sum Of Squares
	* anv[3] - degrees of freedom of ESS
	* anv[4] - degrees of freedom of RSS
	* anv[5] - F Statistics = (R^2/p-1)/((1-R^2)/N-p)
	* anv[6] - P value associated with anv[5] used to reject/accept
	* zero hypothesis
	*/

	if (intercept == 1) {
		if (p == 1) {
			anv[1] = R2[0] * (*bt - N*ym2);
			anv[2] = (1. - R2[0]) * (*bt - N*ym2);
			anv[0] = anv[1] + anv[2];
			anv[3] = (double)p - dfi;
			anv[4] = (double)df;
		}
		else if (p > 1) {
			anv[1] = R2[0] * (*bt - N*ym2);
			anv[2] = (1. - R2[0]) * (*bt - N*ym2);
			anv[0] = anv[1] + anv[2];
			anv[3] = (double)p - dfi;
			anv[4] = (double)df;
			anv[5] = (R2[0] / anv[3]) / ((1. - R2[0]) / anv[4]);
			anv[6] = 1.0 - fcdf(anv[5], p - 1, df);
		}
	}
	else {
		R2[0] = mss / (mss + rss);
		R2[1] = 1.0 - (1.0 - R2[0]) * ((double)N - dfi) / ((double)df);
		if (p == 1) {
			anv[1] = mss;
			anv[2] = rss;
			anv[0] = anv[1] + anv[2];
			anv[3] = (double)p;
			anv[4] = (double)df;
			anv[5] = (R2[0] / anv[3]) / ((1. - R2[0]) / anv[4]);
			anv[6] = 1.0 - fcdf(anv[5], p, df);
		}
		else if (p > 1) {
			anv[1] = mss;
			anv[2] = rss;
			anv[0] = anv[1] + anv[2];
			anv[3] = (double)p;
			anv[4] = (double)df;
			anv[5] = (R2[0] / anv[3]) / ((1. - R2[0]) / anv[4]);
			anv[6] = 1.0 - fcdf(anv[5], p - 1, df);
		}
		else if (p == 0) {
			anv[1] = mss;
			anv[2] = rss;
			anv[0] = anv[1] + anv[2];
			anv[3] = (double)p;
			anv[4] = (double)df;
			anv[5] = NAN;
			anv[6] = NAN;
		}
	}


	free(ipiv);
	free(x);
	free(xt);
	free(xxt);
	free(xxt2);
	free(yxt);
	free(temp1);
	free(temp2);
	free(bt);
	free(seb);
	free(xb);
}

void zerohyp_multi(int N, double *b, int p, double *varcovar, double *tval, double *pval) {
	int df, i;
	/*
	* Zero Hypothesis Test for CLRM
	*/

	df = N - p;

	for (i = 0; i < p; ++i) {
		tval[i] = fabs(b[i]) / sqrt(varcovar[(p + 1)*i]);
		pval[i] = 1.0 - tcdf(fabs(tval[i]), df);
	}

}

void zerohyp_val(reg_object obj, double *tval, double *pval) {
	int df, i;
	/*
	* Zero Hypothesis Test for CLRM
	*/

	df = obj->N - obj->p;
	printf("\n");
	printf("%-25s \n", "Zero Hypothesis Test t values and Probabilities (One-tailed)");
	printf("%-25s%-20s%-20s \n", "Coefficients", "t Value", "Probability");
	if (df > 0) {
		for (i = 0; i < obj->p; ++i) {
			tval[i] = fabs((obj->beta + i)->value) / (obj->beta + i)->stdErr;
			pval[i] = 1.0 - tcdf(fabs(tval[i]), df);
			printf("B%-25d%-20lf%-20g \n", i, tval[i], pval[i]);
		}
	}
	else {
		printf("degrees of freedom cannot be non-positive");
	}
}

void setIntercept(reg_object obj, int intercept) {
	if (intercept == 1 && obj->p > 0) {
		obj->intercept = 1;
	}
	else if (intercept == 0) {
		obj->intercept = 0;
		//obj->p = obj->p - 1;
	}
	else {
		printf("The variable intercept only accepts 0 or 1 values. Additionally, Intercept is 0 when p == 0. \n");
		exit(-1);
	}
}

void setLLSMethod(reg_object obj, char *llsmethod) {
	if (!strcmp(llsmethod, "qr")) {
		strcpy(obj->lls, llsmethod);
	}
	else if (!strcmp(llsmethod, "normal")) {
		strcpy(obj->lls, llsmethod);
		printf("Warning - This method will not calculate the rank of the regression. Use method qr instead. \n");
	}
	else if (!strcmp(llsmethod, "svd")) {
		strcpy(obj->lls, llsmethod);
		printf("Warning - This method will not calculate the rank of the regression. use method qr instead. \n");
	}
	else {
		printf("This function only accepts one of three least square methods - \n");
		printf(" qr, normal and svd. \n");
		exit(-1);
	}
}

void regress(reg_object obj, double *x, double *y, double *res, double *varcovar, double alpha) {
	double *anv2, *b2, *low, *up, *sigma2;
	int p, i, p2, k, dfm;
	double ssr, N2, pi;
	/*
	* obj - Regression object created by reg_init
	*  x - Vector containing (p - 1 ) independent regression variables of length N each
	*  y - Vector containing dependent regression variable of length N
	*  res - Residual vector of length N
	* varcovar - Variance-Covariance Matrix of p regression parameters. Dimension p X p
	* alpha is used to determine confidence interval limits
	* for a given 100*(1-alpha) % confidence interval
	* eg., alpha = 0.05 if you want to determine 95% confidence interval values
	*/
	p = obj->p;
	obj->alpha = alpha;
	anv2 = (double*)malloc(sizeof(double) * 7);
	b2 = p == 0 ? NULL : (double*)malloc(sizeof(double) * p);
	low = (double*)malloc(sizeof(double) * (p + 1));
	up = (double*)malloc(sizeof(double) * (p + 1));
	sigma2 = (double*)malloc(sizeof(double) * 1);

	pi = 3.141592653589793;


	linreg_multi(p, x, y, obj->N, b2, sigma2, varcovar,
		obj->R2, res, alpha, anv2, low, up, &obj->rank, obj->lls, obj->intercept);
	obj->df = obj->N - obj->p;

	obj->sigma = sigma2[0];
	obj->sigma_lower = low[p];
	obj->sigma_upper = up[p];

	obj->TSS = anv2[0];
	obj->ESS = anv2[1];
	obj->RSS = anv2[2];
	obj->df_ESS = (int)anv2[3];
	obj->df_RSS = (int)anv2[4];
	obj->FStat = anv2[5];
	obj->PVal = anv2[6];
	obj->r2 = obj->R2[0];
	obj->r2adj = obj->R2[1];

	for (i = 0; i < p; ++i) {
		(obj->beta + i)->value = b2[i];
		(obj->beta + i)->lower = low[i];
		(obj->beta + i)->upper = up[i];
		(obj->beta + i)->stdErr = sqrt(varcovar[(p + 1)*i]);
	}

	// Log Likelihood

	N2 = (double)(obj->N) / 2.0;

	ssr = 0.0;

	for (i = 0; i < obj->N; ++i) {
		ssr += (res[i] * res[i]);
	}

	obj->loglik = -N2 * log(2 * pi) - N2 * log(ssr / (double)obj->N) - N2;

	k = obj->intercept == 1 ? 1 : 0;
	dfm = obj->intercept == 1 ? p - 1 : p;

	obj->aic = -2.0 * obj->loglik + 2.0 * (double)(k + dfm);
	obj->bic = -2.0 * obj->loglik + log((double)obj->N) * (double)(k + dfm);
	obj->aicc = obj->aic + 2.0 * dfm * ((double)obj->N / ((double)obj->N - dfm - 1.0) - 1.0);

	free(anv2);
	free(b2);
	free(up);
	free(low);
	free(sigma2);
}

void regress_poly(reg_object obj, double *x, double *y, double *res, double *varcovar, double alpha) {
	int polydeg, i, j, sum;
	double *xx;
	/*
	* Performs polynomial regression of the form
	* Y = B0 + B1 * X + B2 * X^2 + ..... + B(p-1) * X ^(p-1)
	*
	* where p is defined as in regular regression
	* p = Number of dependent variables + Independent Variables
	*
	* For p = 3
	* regress_poly calculates coefficients of
	*
	* Y = B0 + B1 * X + B2 * X^2
	*
	*For more information see the functions regress and reg_init
	*/
	if (obj->intercept == 1) {
		polydeg = obj->p - 1;
	}
	else {
		polydeg = obj->p;
	}


	if (polydeg < 1) {
		printf("The Value of p should be greater than or equal to 2");
		exit(1);
	}

	xx = (double*)malloc(sizeof(double) * obj->N * polydeg);
	for (i = 0; i < obj->N; i++)
		xx[i] = x[i];
	sum = obj->N;
	for (i = 1; i < polydeg; ++i) {
		for (j = 0; j < obj->N; ++j) {
			xx[sum + j] = (int)pow((double)x[j], (double)i + 1.);
		}
		sum += obj->N;
	}
	regress(obj, xx, y, res, varcovar, alpha);

	free(xx);
}


void summary(reg_object obj) {
	int i;
	printf("\n");
	printf("%-25s \n", "Regression Summary");
	printf("%-25s%-20s%-20s \n", "Coefficients", "Value", "Standard Error");
	for (i = 0; i < obj->p; ++i) {
		printf("B%-25d%-20lf%-20g \n", i, (obj->beta + i)->value, (obj->beta + i)->stdErr);
	}
	printf("\n");
	printf("Residual Variance = %lf \n", obj->sigma);
	printf("R-Squared = %lf , Adjusted R-Squared = %lf \n", obj->R2[0], obj->R2[1]);

}

void anova(reg_object obj) {
	printf("\n");
	printf("ANOVA : \n");
	printf("%-25s%-20s%-20s%-20s \n", "Source of Variation", "df", "SS", "MSS");
	printf("%-25s%-20d%-20lf%-20lf \n", "Due to Regression",
		obj->df_ESS, obj->ESS, obj->ESS / ((double)obj->df_ESS));
	printf("%-25s%-20d%-20lf%-20lf \n", "Due to Residual",
		obj->df_RSS, obj->RSS, obj->RSS / ((double)obj->df_RSS));
	printf("%-25s%-20d%-20lf \n", "Total",
		obj->df_RSS + obj->df_ESS, obj->TSS);
	printf("\n F Statistics = %g \n", obj->FStat);
	printf("\n P Value (F)  = %g \n", obj->PVal);
}

void anova_list(reg_object *list, int N) {
	int i;
	printf("\n");
	printf("ANOVA Multiple Models : \n");
	printf("%-10s%-25s%-20s%-20s%-20s \n","Model", "Res.Df", "RSS", "Df", "Sum of Sq");
	printf("%-10d%-25d%-20lf \n", 0,
		list[0]->df_RSS, list[0]->RSS);
	for(i = 1; i < N;++i) {
		printf("%-10d%-25d%-20lf%-20d%-20lf \n", i+1,
		list[i]->df_RSS, list[i]->RSS, list[i-1]->df_RSS-list[i]->df_RSS,list[i-1]->RSS-list[i]->RSS);
	}
}

void confint(reg_object obj) {
	int i;
	printf("\n");
	printf("%-10lf%% Confidence Interval For Parameters And Residual Variance : \n", (1.0 - obj->alpha) * 100);
	printf("%-25s%-20s%-20s%-20s \n", "Parameters",
		"Value", "Lower Limit", "Upper Limit");

	for (i = 0; i < obj->p; ++i) {
		printf("B%-25d%-20lf%-20lf%-20lf \n",
			i, (obj->beta + i)->value, (obj->beta + i)->lower, (obj->beta + i)->upper);
	}
	printf("%-25s%-20lf%-20lf%-20lf \n", "Residual Variance",
		obj->sigma, obj->sigma_lower, obj->sigma_upper);
}

double fitted(reg_object obj, double *inp, double *varcovar, double *var) {
	int i, p;
	double oup;
	double *y, *x0, *b, *temp;
	/*
	* The function predicts output value at a given input point (vector of size p-1)
	*
	* inp - ia a  vector of length p-1 values corresponding to each of the p-1
	* dependent variables.
	* The function returns output value corresponding to the input inp vector of length (p-1).
	*
	* varcovar - Is the Aurovariance-Covariance matrix that is obtained from the regression function
	* var - is a vector of length 2 where
	* var[0] - Variance of Mean Prediction
	* var[1] - Variance Of Individual Prediction
	*
	*/
	p = obj->p;
	y = (double*)malloc(sizeof(double) * 1);
	x0 = (double*)malloc(sizeof(double) * p);
	b = (double*)malloc(sizeof(double) * p);
	temp = (double*)malloc(sizeof(double) * p);

	if (obj->intercept == 1) {
		x0[0] = 1.0;

		b[0] = (obj->beta)->value;

		for (i = 1; i < p; ++i) {
			x0[i] = inp[i - 1];
			b[i] = (obj->beta + i)->value;
		}
	}
	else {

		for (i = 0; i < p; ++i) {
			x0[i] = inp[i];
			b[i] = (obj->beta + i)->value;
		}
	}


	// Output
	mmult(x0, b, y, 1, p, 1);
	oup = *y;

	// Variance Of Mean Prediction
	mmult(x0, varcovar, temp, 1, p, p);
	mmult(temp, x0, y, 1, p, 1);
	var[0] = *y;

	// Variance Of Individual prediction
	var[1] = obj->sigma + var[0];

	free(temp);
	free(y);
	free(x0);
	free(b);
	return oup;
}

void free_reg(reg_object object) {
	free(object);

}
/*
void linreg_multi2(int p, double *xi,double *y, int N, double* b) {
double *x,*xt;
double *xxt,*yxt;
double *WORK;
int i,ra,ca,rb,cb;
lapack_int info,lda,infoi,LWORK;
int ipiv[p];
int MONE=-1;
double LWKOPT;


x = (double*) malloc (sizeof(double) * p * N);
xt = (double*) malloc (sizeof(double) * p * N);
xxt = (double*) malloc (sizeof(double) * p * p);
yxt = (double*) malloc (sizeof(double) * p);

// Initalize x matrix
for (i = 0; i < N;++i) {
x[i] = 1.;
}

for (i = N; i < p*N;++i) {
x[i] = xi[i-N];
}

// Calculate x * x'
ra = p;
ca = N;
rb = N;
cb = p;

mtranspose(x,p,N,xt);

cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ra,cb,ca,1.0,x,ca,xt,cb,0.0,xxt,cb);

// Calculate y * x'
ra = 1;
ca = N;
rb = N;
cb = p;

cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ra,cb,ca,1.0,y,ca,xt,cb,0.0,yxt,cb);


// Solve linear equations to find b

// Step 1 Use DGETRF to find the LU factorization

lda = p;

//info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,p,p,xxt,lda,ipiv);

LAPACK_dgetrf(&p,&p,xxt,&lda,ipiv,&info);

//infoi = LAPACKE_dgetri_work(LAPACK_ROW_MAJOR,p,xxt,lda,ipiv,&LWKOPT,MONE);

LAPACK_dgetri(&p,xxt,&lda,ipiv,&LWKOPT,&MONE,&infoi);

LWORK = (int) LWKOPT;

WORK =(double*) malloc (sizeof(double) * LWORK);

//infoi = LAPACKE_dgetri_work(LAPACK_ROW_MAJOR,p,xxt,lda,ipiv,WORK,LWORK);

LAPACK_dgetri(&p,xxt,&lda,ipiv,WORK,&LWORK,&infoi);

// Find b
ra = 1;
ca = p;
rb = p;
cb = p;
cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ra,cb,ca,1.0,yxt,ca,xxt,cb,0.0,b,cb);

free(x);
free(xt);
free(xxt);
free(yxt);
free(WORK);
}
*/

