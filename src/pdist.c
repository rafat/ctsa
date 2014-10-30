/*
 * pdist.c
 *
 *  Created on: Jun 19, 2013
 *      Author: Rafat Hussain
 */

#include "pdist.h"

double normalpdf(double x, double mu, double sigma) {
	double oup,t;
	
	if (sigma < 0.) {
		printf("Standard Deviation (Sigma) must be non-negative");
		exit(1);
	}
	
	t = x - mu;
	t = t * t;
	
	oup = exp(-t / ( 2* sigma *sigma)) / ( sqrt(2 * PIVAL) * sigma);
	
	return oup;
}

double normalcdf(double x, double mu, double sigma) {
	double oup,t;
	
	if (sigma < 0.) {
		printf("Standard Deviation (Sigma) must be non-negative");
		exit(1);
	}
	
	t = (x - mu ) /(sqrt(2.0) * sigma);
	
	oup = 0.5 * (1. + erf(t));
	return oup;
}

double normalinv(double p, double mu, double sigma) {
	double oup,t;
	if (sigma < 0.) {
		printf("Standard Deviation (Sigma) must be non-negative");
		exit(1);
	}
	
	if (p < 0. || p > 1.0) {
		printf("Probablity Values can only take values between 0.0 and 1.0");
		exit(1);
	}
	
	t = sqrt(2.0) * erfinv(2.0 * p - 1.0);
	oup = t * sigma + mu;
	
	return oup;
}

double tpdf(double t, int df) {
	double ft,f1,f2,f3;
	double pi  = 3.1415926535897932384626434;
	
	if (df <= 0) {
		printf("Degree Of Freedom should be a positive real integer");
		exit(1);
	}
	
	f1 = gamma_log((double) (df+1.0)/2.0) - gamma_log((double) (df/2.0));
	f2 = exp(f1);
	f3 = pow(1.0 + (t*t / (df * 1.0)), - (df+1.0)/2.0);
	ft = f2 * f3 ;
	ft = ft / sqrt( (double) pi * df); 
	
	return ft;
	
}

double tcdf(double t, int df) {
	double oup,x,I;
	
	if (df <= 0) {
		printf("Degree Of Freedom should be a positive real integer");
		exit(1);
	}
	
	x = df / (df + t * t);
	
	I = 0.5 * ibeta(x,df / 2.0, 1./2.);
	
	oup = 1. - I;
	
	if ( t < 0.) {
		oup = 1. - oup;
	}
	
	if ( x == 0.) {
		oup = 0.5;
	}
	
	
	return oup;
}

double tinv_appx(double p, int df) {
	double t,x,d,g1,g2,g3,g4;
	/*
	 * Abramowitz and Stegun Approximation 26.7.5
	 */ 
	 
	if (df <= 0) {
		printf("Degree Of Freedom should be a real integer");
		exit(1);
	}
	
	if (p < 0. || p > 1.0) {
		printf("Probablity Values can only take values between 0.0 and 1.0");
		exit(1);
	} 
	
	if ( p == 0.) {
		t = -XINFVAL;
	}
	
	if ( p == 1.) {
		t = XINFVAL;
	}
	
	if (p > 0. && p < 1.0) {
		x = normalinv(p,0.0,1.0);
		d = (double) df;
		g1 = (pow(x,3.0) + x ) / (4. * d);
		g2 = (5. * pow(x,5.0) + 16. * pow(x,3.0) + 3. * x) / (96. * pow(d, 2.0)); 
		g3 = (3. * pow(x,7.0) + 19. * pow(x,5.0) + 17. * pow(x,3.0) - 15. * x)\
		/(384. * pow(d, 3.0)); 
		g4 = (79. * pow(x,9.0) + 776. * pow(x,7.0) + 1482. * pow(x,5.0) - 1920. *\
		 pow(x,3.0) - 945. * x)/(92160. * pow(d, 4.0)); 
		 t = x + g1 + g2 + g3 + g4;
	}
	
	return t;
}

double tinv(double p, int df) {
	double t,q,pmin,x,y,d2,sign;
	
	if (df <= 0) {
		printf("Degree Of Freedom should be a positive real integer");
		exit(1);
	}
	
	if (p < 0. || p > 1.0) {
		printf("Probablity Values can only take values between 0.0 and 1.0");
		exit(1);
	}
	
	q = 1. - p;
	
	if ( p <= q ) {
		pmin = 2.0 * p;
	} else {
		pmin = 2.0 * q;
	}
	
	// XINFVAL value is set to 1.79e+308. It can be changed in erfunc.h file.
	
	if ( p == 0.) {
		t = -XINFVAL;
	}
	
	if ( p == 1.) {
		t = XINFVAL;
	}
	
	
	if ( p > 0. && p < 1. && df <= 1000 ) {
		if (p < 0.5) {
			sign = - 1.0;
		} else {
			sign = 1.0;
		}
		
		d2 = ((double) df ) / 2.0;
		x = betainv(pmin,d2,0.5);
		y = 1. - x;
		
		t = sqrt( (double) df * y / x);
	}
	
	if ( p > 0. && p < 1. && df > 1000 ) {
		t = tinv_appx(p,df);
	}
	
	return t;
}

double fpdf(double x, int num,int den) {
	double oup,k1,k2,y,z;
	
	if (num <= 0 || den <= 0) {
		printf("Degrees Of Freedom should be real integers");
		exit(1);
	}
	
	if (x < 0.) {
		printf("The range of distribution only takes non-negative and real values");
		exit(1);
	}
	
	k1 = (double) num;
	k2 = (double) den;
	
	if ( k1 * x > k2 ) {
		y = (k2 * k1) / ((k2 + k1 * x) * (k2 + k1 * x)) ;
		oup = y * betapdf(k2 / (k2 + k1 * x), k2 / 2.,k1 / 2.);
	} else {
		z = k2 + k1 * x;
		y = (z * k1 - x * k1 * k1) / (pow(z,2.0));
		oup = y * betapdf(k1 * x / (k2 + k1 * x), k1 / 2.,k2 / 2.);
	}
	
	return oup;
}

double fcdf(double x, int num,int den) {
	double oup,k1,k2,y;
	
	if (num <= 0 || den <= 0) {
		printf("Degrees Of Freedom should be positive real integers");
		exit(1);
	}
	
	if (x < 0.) {
		printf("The range of distribution only takes non-negative and real values");
		exit(1);
	}
	
	k1 = (double) num;
	k2 = (double) den;
	
	if (k1 * x < k2) {
		y = k1 * x / (k2 + k1 * x);
		oup = ibeta(y,k1/2.,k2/2.);
	} else {
		y = k2 / (k2 + k1 * x);
		oup = ibetac(y,k2/2.,k1/2.);
	}
	
	return oup;
}

double finv(double p, int num,int den) {
	double x,k1,k2,a,b;
	
	if (p < 0. || p > 1.0) {
		printf("Probablity Values can only take values between 0.0 and 1.0");
		exit(1);
	}
	
	if (num <= 0 || den <= 0) {
		printf("Degrees Of Freedom should be positive real integers");
		exit(1);
	}
	
	k1 = (double) num;
	k2 = (double) den;
	
	a = betainv(p,k1/2.,k2/2.);
	b = 1. - a;
	
	x = (k2 * a) / ( k1 * b);
	
	if ( p == 1.) {
		x = XINFVAL;
	}
	
	return x;
}

double gammapdf(double x, double k, double th) {
	double oup,t;
	/*
	 * Gamma pdf = x ^ (k-1) * exp( - k/th) / (gamma(k) * theta ^k)
	 * Define t = x /th
	 * pdf = exp ((k-1) *log(t) - t - gamma_log(k)) /th
	 */ 
	
	if (k <= 0. || th <= 0.) {
		printf("Gamma Distribution parameters must be positive and real.");
		exit(1);
	}
	
	t = x / th;
	
	if (t > 0.) {
		oup = exp ((k-1) *log(t) - t - gamma_log(k)) /th;
	} else if (t < 0.) {
		oup = 0.0;
	} else if (t == 0.0 && k == 1) {
		oup = 1.0 / th;
	} else {
		printf("The Output is Complex and cannot be determined");
		exit(1);
	}
	
	return oup;
}

double gammacdf(double x, double k, double th) {
	double oup,t;
	
	if (k <= 0. || th <= 0.) {
		printf("Gamma Distribution parameters must be positive and real.");
		exit(1);
	}
	
	if (x <= 0.) {
		x = 0.;
	}
	
	t = x/th;
	
	oup = pgamma(t,k);
	
	return oup;
}

double gammainv(double p, double k, double th) {
	double oup,x,xi,del,mu,sigma2,cvar;
	int m;
	if (p < 0. || p > 1.0) {
		printf("Probablity Values can only take values between 0.0 and 1.0");
		exit(1);
	}
	
	if (k <= 0. || th <= 0.) {
		printf("Gamma Distribution parameters must be positive and real.");
		exit(1);
	}
	
	if (p == 0.) {
		oup = 0.0;
	} 
	
	if ( p == 1.) {
		oup = XINFVAL;
	}
	
	if ( p > 0.0 && p < 1.0) {
		// Use Newton's Method to find cdf inverse of gamma distribution
		
		// Initialization is based on lognormal distribution with
		// mean and variance k
		
		sigma2 = log(1+k) - log(k);
		mu = log(k) - 0.5 * sigma2;
		xi = exp(mu - sqrt(2. * sigma2) * erfcinv(2*p));
		
		del = 1.0;
		
		cvar = eps(XNINFVAL);
		m = 0;
		
		
		//Iterations
		
		while ( ( fabs(del) > cvar * xi) &&  (m < 1000)) {
			del = (gammacdf(xi,k,1.0) - p) / r8_max(gammapdf(xi,k,1.0), XNINFVAL);
			//printf("%g %g \n",xi*cvar,del);
			m++;
			x = xi - del;
			
			if (x <= 0.) {
				x = xi / 10.;
				del = xi - x;
			}
			xi = x;
		}
		oup = xi * th;
	}
	
	return oup;
}

double chipdf(double x, int df) {
	double oup,k;
	if (df <= 0 ) {
		printf("Degree of freedom is a real positive integer");
		exit(1);
	}
	
	if (x < 0.) {
		printf("The range of distribution only takes non-negative and real values");
		exit(1);
	}
	
	k = (double) df;
	oup = gammapdf(x,k/2,2.);
	
	return oup;
}

double chicdf(double x, int df) {
	double oup,a;
	if (df <= 0 ) {
		printf("Degree of freedom is a real positive integer");
		exit(1);
	}
	
	if (x < 0.) {
		printf("The range of distribution only takes non-negative and real values");
		exit(1);
	}
	a = (double) df;
	oup = pgamma(x/2,a/2);
	return oup;
}

double chiinv(double p, int df) {
	double x,k;
	
	if (df <= 0 ) {
		printf("Degree of freedom is a real positive integer");
		exit(1);
	}
	
	if (p < 0. || p > 1.0) {
		printf("Probablity Values can only take values between 0.0 and 1.0");
		exit(1);
	}
	
	k = (double) df;
	
	x = gammainv(p,k/2.,2.0);
	
	return x;
}
