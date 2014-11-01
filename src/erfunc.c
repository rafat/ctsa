/*
 *
 * This program is modified C translation of ACM TOMS 708 erf and erfc
 * functions. This translation was done by R Project
 * http://svn.r-project.org/R/trunk/src/nmath/
 * The original, true, correct version of ACM TOMS Algorithm 708 is
 * available through ACM: http://www.acm.org/pubs/calgo or
 * NETLIB: http://www.netlib.org/toms/index.html.
 *
 *    References
 *
 *  1.  Armido Didonato, Alfred Morris, Jr,
 *	Algorithm 708: Significant Digit Computation of the Incomplete Beta Function Ratios,
 *	ACM Transactions on Mathematical Software,
 *  Volume 18, Number 3, pages 360-373, 1992.
 */

#include "erfunc.h"

static double exparg(int l)
{
/* -------------------------------------------------------------------- */
/*     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
 *     EXP(W) CAN BE COMPUTED.  ==>  exparg(0) = 709.7827  nowadays. */

/*     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
 *     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
 *       ==> exparg(1) = -708.3964   nowadays. */

/*     Note... only an approximate value for exparg(L) is needed. */
/* -------------------------------------------------------------------- */

    double lnb = .69314718055995;
    int m = (l == 0) ? DBL_MAX_EXP : DBL_MIN_EXP - 1;

    return m * lnb * .99999;
} /* exparg */

double erf__(double x)
{
/* -----------------------------------------------------------------------
 *             EVALUATION OF THE REAL ERROR FUNCTION
 * ----------------------------------------------------------------------- */

    /* Initialized data */

     double c = .564189583547756;
     double a[5] = { 7.7105849500132e-5,-.00133733772997339,
	    .0323076579225834,.0479137145607681,.128379167095513 };
     double b[3] = { .00301048631703895,.0538971687740286,
	    .375795757275549 };
     double p[8] = { -1.36864857382717e-7,.564195517478974,
	    7.21175825088309,43.1622272220567,152.98928504694,
	    339.320816734344,451.918953711873,300.459261020162 };
     double q[8] = { 1.,12.7827273196294,77.0001529352295,
	    277.585444743988,638.980264465631,931.35409485061,
	    790.950925327898,300.459260956983 };
     double r[5] = { 2.10144126479064,26.2370141675169,
	    21.3688200555087,4.6580782871847,.282094791773523 };
     double s[4] = { 94.153775055546,187.11481179959,
	    99.0191814623914,18.0124575948747 };

    /* System generated locals */
    double ret_val;

    /* Local variables */
    double t, x2, ax, bot, top;

    ax = fabs(x);
    if (ax <= 0.5) {
	t = x * x;
	top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0;
	bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;

	return x * (top / bot);
    }
    /* else: ax > 0.5 */

    if (ax <= 4.) { /*  ax in (0.5, 4] */

	top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
		+ p[5]) * ax + p[6]) * ax + p[7];
	bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
		+ q[5]) * ax + q[6]) * ax + q[7];
	ret_val = 0.5 - exp(-x * x) * top / bot + 0.5;
	if (x < 0.0) {
	    ret_val = -ret_val;
	}
	return ret_val;
    }

    /* else: ax > 4 */

    if (ax >= 5.8) {
	return x > 0 ? 1 : -1;
    }
    x2 = x * x;
    t = 1.0 / x2;
    top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
    bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
    t = (c - top / (x2 * bot)) / ax;
    ret_val = 0.5 - exp(-x2) * t + 0.5;
    if (x < 0.0) {
	ret_val = -ret_val;
    }
    return ret_val;

} /* erf */

double erfc1__(int ind, double x)
{
/* ----------------------------------------------------------------------- */
/*         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION */

/*          ERFC1(IND,X) = ERFC(X)            IF IND = 0 */
/*          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE */
/* ----------------------------------------------------------------------- */

    /* Initialized data */

    double c = .564189583547756;
    double a[5] = { 7.7105849500132e-5,-.00133733772997339,
	    .0323076579225834,.0479137145607681,.128379167095513 };
    double b[3] = { .00301048631703895,.0538971687740286,
	    .375795757275549 };
    double p[8] = { -1.36864857382717e-7,.564195517478974,
	    7.21175825088309,43.1622272220567,152.98928504694,
	    339.320816734344,451.918953711873,300.459261020162 };
    double q[8] = { 1.,12.7827273196294,77.0001529352295,
	    277.585444743988,638.980264465631,931.35409485061,
	    790.950925327898,300.459260956983 };
    double r[5] = { 2.10144126479064,26.2370141675169,
	    21.3688200555087,4.6580782871847,.282094791773523 };
    double s[4] = { 94.153775055546,187.11481179959,
	    99.0191814623914,18.0124575948747 };

    double ret_val;
    double e, t, w, bot, top;

    double ax = fabs(x);
    //				|X| <= 0.5 */
    if (ax <= 0.5) {
	double t = x * x,
	    top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0,
	    bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;
	ret_val = 0.5 - x * (top / bot) + 0.5;
	if (ind != 0) {
	    ret_val = exp(t) * ret_val;
	}
	return ret_val;
    }
    // else (L10:):		0.5 < |X| <= 4
    if (ax <= 4.0) {
	top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax
		+ p[5]) * ax + p[6]) * ax + p[7];
	bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax
		+ q[5]) * ax + q[6]) * ax + q[7];
	ret_val = top / bot;

    } else { //			|X| > 4
	// L20:
	if (x <= -5.6) {
	    // L50:            	LIMIT VALUE FOR "LARGE" NEGATIVE X
	    ret_val = 2.0;
	    if (ind != 0) {
		ret_val = exp(x * x) * 2.0;
	    }
	    return ret_val;
	}
	if (ind == 0 && (x > 100.0 || x * x > -exparg(1))) {
	    // LIMIT VALUE FOR LARGE POSITIVE X   WHEN IND = 0
	    // L60:
	    return 0.0;
	}

	// L30:
	t = 1. / (x * x);
	top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
	bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
	ret_val = (c - t * top / bot) / ax;
    }

    // L40:                 FINAL ASSEMBLY
    if (ind != 0) {
	if (x < 0.0)
	    ret_val = exp(x * x) * 2.0 - ret_val;
    } else {
	// L41:  ind == 0 :
	w = x * x;
	t = w;
	e = w - t;
	ret_val = (0.5 - e + 0.5) * exp(-t) * ret_val;
	if (x < 0.0)
	    ret_val = 2.0 - ret_val;
    }
    return ret_val;

} /* erfc1 */


double erf(double x) {
	return erf__(x);
}

double erfc(double x) {
	return erfc1__(0,x);
}

double erfcx(double x) {
	// returns EXP(X*X)*ERFC(X)
	return erfc1__(2,x);

}

double erfinv(double x) {
	double a,t,oup,num,den,xinf,sign,temp,pi;
	/*
	 * This implementation is based on 
	 * J. M. Blair, C. A. Edwards, and J. H. Johnson, "Rational
	Chebyshev Approximations for the Inverse of the Error Function",
	Mathematics of Computation, 30 (1976) 827-830 
	* 
	* Tables from www.jstor.org's copy of the paper are used in this
	 * implementation.
	 */ 
	 if ( x > 0.0) {
		 sign = 1.0; 
	} else {
		 sign = -1.0;	 
	}
	 a = fabs(x);
	 xinf = 1.79e308;
	 pi  = 3.1415926535897932384626434;
	 
	 if (a > 1.) {
		 printf("The Input takes values between -1 and +1");
		 exit(1);
	 }
	 
	 if (a <= 0.75) {
		 //Coefficients from table 13
		 double p[5] = {4.62680202125696,-16.6805947126248,17.6230176190819,\
			 -5.4309342907266,0.236997019142};
		 double q[5] = {4.26606447606664,-17.5930726990431,22.7331464544494,\
			 -9.9016863476727,1.000000000000};
		 
		 t = x*x - 0.75*0.75;
		 num = p[0] + t * ( p[1] + t * ( p[2] + t * ( p[3] + t * p[4]))); 
		 den = q[0] + t * ( q[1] + t * ( q[2] + t * ( q[3] + t * q[4]))); 
		 oup = x * num / den;
		 
	 } else if ( 0.75 < a && a <= 0.9375) {
		 //Coefficients from table 33
		 double p[6] = {-0.041199817067782,0.643729854003468,-3.28902674093993,\
			6.24518431579026,-3.65953596495117,0.30260114943200};
		 double q[6] = {-0.029324540620124,0.501148500527886,-2.90144687299145,\
			 6.65393051963183,-5.40640580412825,1.0000000000000};
		 
		 t = x*x - 0.9375*0.9375;
		 num = p[0] + t * ( p[1] + t * ( p[2] + t * ( p[3] + t * (p[4] + t * p[5])))); 
		 den = q[0] + t * ( q[1] + t * ( q[2] + t * ( q[3] + t * (q[4] + t * q[5])))); 
		 oup = x * num / den;	 
	 } else if ( 0.9375 < a && a  < 1.) {
		 //Coefficients from table 53
		 /*double p[7] = {0.00479447793346489,0.17121766901388401,1.26186520719596240,\
			 2.1339060556433715,0.2838737490019514,0.2353035904898929,\
			 -0.0400763898644416};
		 double q[6] = {0.00479437046609729,0.17124542554915192,1.27422170654959402,\
			 2.4161883799488509,1.4957337483050566,1.0000000000000000};
			 
		 t=1.0/sqrt(-log(1.0-a));
		 num = p[0] + t * ( p[1] + t * ( p[2] + t * ( p[3] + t * ( p[4] + t * ( p[5] + t * ( p[6] + t * p[7]))))));	 
		 den = q[0] + t * ( q[1] + t * ( q[2] + t * ( q[3] + t * (q[4] + t * ( q[5] + t * q[6]))))); 
		 oup = sign * num / (den * t);	
		 * */
		 
		//Coefficients from table 50 as these are more accurate
		double p[6] = {.1550470003116,1.382719649631,.690969348887, \
			-1.128081391617, .680544246825,-.16444156791};
		double q[3] = {.155024849822,1.385228141995,1.000000000000}; 
		
		t=1.0/sqrt(-log(1.0-a));
		oup = sign*(p[0]/t+p[1]+t*(p[2]+t*(p[3]+t*(p[4]+t*p[5]))))/
          (q[0]+t*(q[1]+t*(q[2])));
	 } else {
		 oup = sign * xinf;
	 } 
	 
	 // Apply Newton's correction and Halley's step
	 
	 if (a < 1.) {
		 temp = (erf(oup) - x) / ( 2 / sqrt(pi) * exp(- pow(oup,2.0)));
		 oup -= temp / ( 1 + temp*oup);
	 }
	 
	return oup;
}

double erfcinv(double x) {
	double oup,xinf,pi,t;
	
	/*
	 * This implementation is based on 
	 * J. M. Blair, C. A. Edwards, and J. H. Johnson, "Rational
	Chebyshev Approximations for the Inverse of the Error Function",
	Mathematics of Computation, 30 (1976) 827-830 
	* 
	* Tables from www.jstor.org's copy of the paper are used in this
	 * implementation.
	 */ 
	
	xinf = 1.79e308;
	pi  = 3.1415926535897932384626434;
	
	if (x > 2. || x < 0.) {
		 printf("The Input takes values between 0.0 and 2.0");
		 exit(1);
	 }
	 
	if( x >= 0.0625 && x <2.0) {
		oup = erfinv( 1.0 - x );
	} else if (x < 0.0625 && x >= 1.0e-100) {
		//Coefficients from table 50 
		double p[6] = {0.1550470003116,1.382719649631,0.690969348887, \
			-1.128081391617, 0.680544246825,-0.16444156791};
		double q[3] = {0.155024849822,1.385228141995,1.000000000000};
		
		t=1.0/sqrt(-log(x));
		oup = (p[0]/t+p[1]+t*(p[2]+t*(p[3]+t*(p[4]+t*p[5]))))/
          (q[0]+t*(q[1]+t*(q[2])));
		
	} else if( x < 1.0e-100 && x > 1.0e-1000 ) {
		//Coefficients from table 70
		double p[4]={0.00980456202915,0.363667889171,0.97302949837,-0.5374947401};
		double q[3]={0.00980451277802,0.363699971544,1.000000000000}; 
		t = 1.0/sqrt(-log(x));
		oup = (p[0]/t+p[1]+t*(p[2]+t*p[3]))/(q[0]+t*(q[1]+t*(q[2])));
		
	} else if (x <= 1.0e-1000) {
		 oup = xinf;
	} else if (x == 2.) {
		 oup = -xinf;
	}
	
	return oup;
	
}
