/*
 * dist.c

 *
 *  Created on: Jun 11, 2013
 *      Author: Rafat Hussain
 */

#include "dist.h"

double fix(double x) {
	// Rounds to the integer nearest to zero 
	if (x >= 0.) {
		return floor(x);
	} else {
		return ceil(x);
	}
}

double normf(double x) {
	if ( x < 0.0) {
		x = - x;
	}
	return x;
}

double eps(double x) {
	double y,t,epsv;
	int t2;
	
	epsv = pow(2.0,-52.0);
	if (x < 0.) {
		x = - x;
	}
	
	if ( x > 0. && x < 1.0) {
		x /= 2.;
	}
	t2 = (int) ( log(x) / log ((double) 2.0));
	t = pow(2.0,(double) t2);
	
	y = epsv * t;
	return y;
}

double log1p(double x) {
	/*
	 * This log(1+x) algorithm is courtesy of Charles kerney and John D Cook
	 * http://www.johndcook.com/blog/2012/07/25/trick-for-computing-log1x/
	 * and 
	 * http://geographiclib.sourceforge.net/html/Math_8hpp_source.html#l00241
	 * 
	 */ 
	volatile double y,z;
	y = 1 + x;
	z = y - 1;
	return z == 0 ? x : x * log(y) / z;
}

double recfact(double N) {
	double i;

	i = (int) N;

	if (i == 1) {
		return 1;
	} else {
		return N * recfact(N-1);
	}

}

double factorial(int N) {
	double i,ans;
	i = (double) N;
	ans = recfact(i);
	return ans;
}

double gamma(double x) {
	/*
	 * This C program code is based on  W J Cody's fortran code.
	 * http://www.netlib.org/specfun/gamma
	 * 
	 * References:
   "An Overview of Software Development for Special Functions",
	W. J. Cody, Lecture Notes in Mathematics, 506,
	Numerical Analysis Dundee, 1975, G. A. Watson (ed.),
	Springer Verlag, Berlin, 1976.

   Computer Approximations, Hart, Et. Al., Wiley and sons, New York, 1968.
   */ 
   
	// numerator and denominator coefficients for 1 <= x <= 2
	
	double y,oup,fact,sum,y2,yi,z,nsum,dsum;
	int swi,n,i;
	
	double spi = 0.9189385332046727417803297;
    double pi  = 3.1415926535897932384626434;
	double xmax = 171.624e+0;
	double xinf = 1.79e308;
	double eps = 2.22e-16;
	double xninf = 1.79e-308;
	
	double num[8] = { -1.71618513886549492533811e+0,
                        2.47656508055759199108314e+1,
                       -3.79804256470945635097577e+2,
                        6.29331155312818442661052e+2,
                        8.66966202790413211295064e+2,
                       -3.14512729688483675254357e+4,
                       -3.61444134186911729807069e+4,
                        6.64561438202405440627855e+4 };
                       
    double den[8] = { -3.08402300119738975254353e+1,
                        3.15350626979604161529144e+2,
                       -1.01515636749021914166146e+3,
                       -3.10777167157231109440444e+3,
                        2.25381184209801510330112e+4,
                        4.75584627752788110767815e+3,
                       -1.34659959864969306392456e+5,
                       -1.15132259675553483497211e+5 }; 
                       
    // Coefficients for Hart's Minimax approximation x >= 12       
    
    
	double c[7] = { -1.910444077728e-03,
                        8.4171387781295e-04,
                       -5.952379913043012e-04,
                        7.93650793500350248e-04,
                       -2.777777777777681622553e-03,
                        8.333333333333333331554247e-02,
                        5.7083835261e-03 };            
    
    y = x;             
    swi = 0; 
    fact = 1.0;
    n = 0;      
    
    
    if ( y < 0.) {
		// Negative x
		y = -x;
		yi = fix(y);
		oup = y - yi;
		
		if (oup != 0.0) {
			if (yi != fix(yi * .5) * 2.) {
				swi = 1;
			}
			fact = -pi / sin(pi * oup);
			y += 1.;
		} else {
			return xinf;
		}
	}      
	
	if (y < eps) {
		if (y >= xninf) {
			oup = 1.0/y;
		} else {
			return xinf;
		}
		
	} else if (y < 12.) {
		yi = y;
		if ( y < 1.) {
			z = y;
			y += 1.;
		} else {
			n = ( int ) y - 1;
            y -= ( double ) n;
            z = y - 1.0;
		}
		nsum = 0.;
		dsum = 1.;
		for (i = 0; i < 8; ++i) {
			nsum = (nsum + num[i]) * z;
			dsum = dsum * z + den[i];
		}
		oup = nsum / dsum + 1.;
		
		if (yi < y) {
	   
			oup /= yi;
		} else if (yi > y) {
	    
			for (i = 0; i < n; ++i) {
				oup *= y;
				y += 1.;
			}
		
		}      
	
	} else {
		if (y <= xmax) {
			y2 = y * y;
			sum = c[6];
			for (i = 0; i < 6; ++i) {
				sum = sum / y2 + c[i];
			}
			sum = sum / y - y + spi;
			sum += (y - .5) * log(y);
			oup = exp(sum);
		} else {
			return(xinf);
		}
	}
	
	if (swi) {
		oup = -oup;
	}
    if (fact != 1.) {
		oup = fact / oup; 
	}       
	
	return oup;
}

double gamma_log(double x) {
	/*
	 * This C program code is based on  W J Cody's fortran code.
	 * http://www.netlib.org/specfun/gamma
	 * 
	 *      References:                                                                                                                                  *
      [1] "Chebyshev Approximations for the Natural Logarithm of the Gamma    
           Function", W. J. Cody and K. E. Hillstrom, Math. Comp. 21, 1967,   
           pp. 198-203.                                                       
                                                                              
      [2] K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May 1969. 
                                                                              
      [3] Hart, Et. Al., Computer Approximations, Wiley and sons, New York,   
          1968.                                                               
   */ 
   
	// numerator and denominator coefficients for 1 <= x <= 2
	
	double y,oup,fact,y2,nsum,dsum;
	int i;
	
	double spi = 0.9189385332046727417803297;
    double pi  = 3.1415926535897932384626434;
	double xmax = 171.624e+0;
	double xinf = 1.79e308;
	double eps = 2.22e-16;
	double xinf4 = 1.1567e+077;
	double pnt = 0.6796875e+0;
	
	double d1 = -5.772156649015328605195174e-1;
	
	// Coefficients for 0.5 < x <= 1.5
	
	double p1[8] = { 4.945235359296727046734888e+0,
                        2.018112620856775083915565e+2,
                        2.290838373831346393026739e+3,
                        1.131967205903380828685045e+4,
                        2.855724635671635335736389e+4,
                        3.848496228443793359990269e+4,
                        2.637748787624195437963534e+4,
                        7.225813979700288197698961e+3 };
                        
    double q1[8] = { 6.748212550303777196073036e+1,
                        1.113332393857199323513008e+3,
                        7.738757056935398733233834e+3,
                        2.763987074403340708898585e+4,
                        5.499310206226157329794414e+4,
                        6.161122180066002127833352e+4,
                        3.635127591501940507276287e+4,
                        8.785536302431013170870835e+3 };            
                        
	// Coefficients for 1.5 < x <= 4.0
	
	double d2 = 4.227843350984671393993777e-1;
	
	double p2[8] = { 4.974607845568932035012064e+0,
                        5.424138599891070494101986e+2,
                        1.550693864978364947665077e+4,
                        1.847932904445632425417223e+5,
                        1.088204769468828767498470e+6,
                        3.338152967987029735917223e+6,
                        5.106661678927352456275255e+6,
                        3.074109054850539556250927e+6 };
                        
    double q2[8] = { 1.830328399370592604055942e+2,
                        7.765049321445005871323047e+3,
                        1.331903827966074194402448e+5,
                        1.136705821321969608938755e+6,
                        5.267964117437946917577538e+6,
                        1.346701454311101692290052e+7,
                        1.782736530353274213975932e+7,
                        9.533095591844353613395747e+6 };                 
    
                        
	// Coefficients for 4.0 < x <= 12.0
	
	double d4 = 1.791759469228055000094023e+0;

	double p4[8] = { 1.474502166059939948905062e+04,
                        2.426813369486704502836312e+06,
                        1.214755574045093227939592e+08,
                        2.663432449630976949898078e+09,
                        2.940378956634553899906876e+10,
                        1.702665737765398868392998e+11,
                        4.926125793377430887588120e+11,
                        5.606251856223951465078242e+11 };
	
	double q4[8] = { 2.690530175870899333379843e+03,
                        6.393885654300092398984238e+05,
                        4.135599930241388052042842e+07,
                        1.120872109616147941376570e+09,
                        1.488613728678813811542398e+10,
                        1.016803586272438228077304e+11,
                        3.417476345507377132798597e+11,
                        4.463158187419713286462081e+11 };
                       
    // Coefficients for Hart's Minimax approximation x >= 12       
    
    
	double c[7] = { -1.910444077728e-03,
                        8.4171387781295e-04,
                       -5.952379913043012e-04,
                        7.93650793500350248e-04,
                       -2.777777777777681622553e-03,
                        8.333333333333333331554247e-02,
                        5.7083835261e-03 };     
	
	oup = x;
	
	if (x < 0) {
		printf("Input must be Positive and Real");
		exit(1);
	}
	
	if (x > xinf) {
		printf("Overflow");
		exit(1);
	}
	
	if (x <= eps) {
		oup = - log(x);
	}
	
	if ( x > eps && x <= 0.5) {
		y = x;
		nsum = 0.0;
		dsum = 1.0;
		
		for (i = 0; i < 8; ++i) {
			nsum = nsum * y + p1[i];
			dsum = dsum * y + q1[i];
		}
		
		oup = -log(y) + (y * (d1 + y * (nsum / dsum))); 
		
	}
	
	if ( x > 0.5 && x <= pnt) {
		y = (x - 0.5) - 0.5;
		nsum = 0.0;
		dsum = 1.0;
		
		for (i = 0; i < 8; ++i) {
			nsum = nsum * y + p2[i];
			dsum = dsum * y + q2[i];
		}
		
		oup = -log(x) + (y * (d2 + y * (nsum / dsum))); 
		
	}
	
	if ( x > pnt && x <= 1.5) {
		y = (x - 0.5) - 0.5;;
		nsum = 0.0;
		dsum = 1.0;
		
		for (i = 0; i < 8; ++i) {
			nsum = nsum * y + p1[i];
			dsum = dsum * y + q1[i];
		}
		
		oup = (y * (d1 + y * (nsum / dsum))); 
		
	}
	
	if ( x > 1.5 && x <= 4.0) {
		y = x - 2.;
		nsum = 0.0;
		dsum = 1.0;
		
		for (i = 0; i < 8; ++i) {
			nsum = nsum * y + p2[i];
			dsum = dsum * y + q2[i];
		}
		
		oup = (y * (d2 + y * (nsum / dsum))); 
		
	}
	
	if ( x > 4.0 && x <= 12.0) {
		y = x - 4.;
		nsum = 0.0;
		dsum = -1.0;
		
		for (i = 0; i < 8; ++i) {
			nsum = nsum * y + p4[i];
			dsum = dsum * y + q4[i];
		}
		
		oup = d4 + y * (nsum / dsum); 
		
	}
	
	if (x > 12.0) {
		oup = 0.0;
		y = x;
		if (y <= xinf4) {
			oup = c[6];
            y2 = y * y;
			for ( i = 0; i < 6; i++ ) {
				oup = oup / y2 + c[i];
			}
		}
		oup /= y;
		fact = log(y);
		oup += spi - 0.5 * fact + y * (fact - 1.0);
		
	}
	
	return oup;
                        

}

double beta(double a, double b) {
	double oup;
	if (a < 0. || b < 0.) {
		printf(" The Two Inputs should be nonnegative and real");
		exit(1);
	}
	oup = exp(gamma_log(a) + gamma_log(b) - gamma_log(a+b)); 
	return oup;
}

double beta_log(double a, double b) {
	double oup;
	if (a < 0. || b < 0.) {
		printf(" The Two Inputs should be nonnegative and real");
		exit(1);
	}
	oup = gamma_log(a) + gamma_log(b) - gamma_log(a+b); 
	return oup;
}


double pgamma(double x, double a) {
	double oup,ct,d,s,la;
	double a0,a1,b0,b1,g,g0,delta,fact;
	// Negative a
	if (a < 0.) {
		printf("a must be non-negative");
		exit(1);
	}
	
	// a -> inf
	la = 2.0e20;
	if (a > la) {
		x = la - 1.0/3.0 + sqrt(la/a)*(x - (a-1.0/3.0));
		a = la;
		if ( x < 0.0) {
			x = 0.0;
		}
	}
	
	// x < a + 1
	
	if ( a != 0 && x != 0 && x < a + 1) {
		ct = a;
		d = 1.;
		s = 1.;
		
		while (normf(d) >= 100*eps(s)) {
			ct += 1.0;
			d = x*d/ct;
			s += d;
		}
		oup = s * exp( -x + a * log(x) - gamma_log(a + 1));
		if ( x > 0 && oup > 1.) {
			oup = 1.0;
		}
	}
	
	if ( a != 0 && x >= a + 1) {
		a0 = 1.0;
		a1 = x;
		b0 = 0.0;
		b1 = 1.0;
		fact = 1.0/a1;
		ct = 1.0;
		g = b1 * fact;
		g0 = b0;
		delta = normf(g - g0);
		
		while (delta >= 100*eps(normf(g))) {
			g0 = g;
			d = ct - a;
			a0 = (a1 + a0*d)*fact;
			b0 = (b1 + b0*d)*fact;
			s = ct*fact;
			a1 = x * a0 + s * a1;
			b1 = x * b0 + s * b1;
			fact = 1.0/a1;
			ct += 1.0;
			g = b1 * fact;
			delta = normf(g - g0);
			
		}
		oup = 1.0 - g * exp( -x + a * log(x) - gamma_log(a));
		
	}
	
	if ( x == 0) {
		oup = 0.0;
	}
	
	if ( a == 0) {
		oup = 1.0;
	}
	
	
	return oup;
}

double qgamma(double x, double a) {
	double oup;
	oup = pgamma(x,a);
	return 1.0 - oup;
	
}

double bfrac(double x,double a, double b, int* ctr) {
	double d2m,d2m1,a0,b0,a1,b1,am,bm,g,g0,delta,a1m,a2m,del;
	int m;
	
	/* This continued fraction computation is based on
	 * equation 26.5.8 Handbook of Mathematical Functions.
	 * by Abramowitz and Stegun. 
	 * If convergence is not achieved in 1000 iterations then
	 * the function will still return values so it is recommended
	 * that calling function 
	 * 
	 * y = bfrac(x,a,b,&controlint)
	 * 
	 * should check controlint and discard the result y
	 * if controlint = 1000.
	 * 
	 */ 
	 m = 1;
	 am = bm = g = 1.;
	 a0 = b0 = a1 = b1 = g0 = 0.;
	 del = 1. - (a+b) * x / (a+1);
	 delta = fabs(g - g0);
	 
	 while ( (delta > 1000*eps(g) ) && (m < 1000) ) {
		 a1m = (double) a + m;
		 a2m = (double) a1m + m;
		 g0 = g;
		 d2m = m * (b - m) * x / (a2m * (a2m - 1.));
		 d2m1 = - a1m * (a1m + b) * x / (a2m * (a2m + 1.));
		 
		 
		 a0 = g + d2m * am;
		 b0 = del + d2m * bm;
		 a1 = a0 + d2m1 * g;
		 b1 = b0 + d2m1 * del;
		 if (m == 1) {
			 del = 1.0;
		 }
		 
		 am = a0 / b1;
		 bm = b0 / b1;
		 g = a1/b1;
		 
		 m++;
		 delta = fabs(g - g0);
	 }
	 
	 *ctr = m;
	 return g;
}

double ibeta_appx(double x, double a , double b) {
	double oup,abx,c2,w1,w2,num,den;
	/*
	 * Method 2. uses approximation and is based on equations 26.5.20 
	 * (pgamma/qgamma) and 26.5.21 (erf)
	 */ 
	 
	 if (x < 0. || x > 1.) {
		printf("x only accepts real values between 0.0 and 1.0");
		exit(1);
	}
	if (a < 0. || b < 0.) {
		printf(" The Two Inputs should be nonnegative and real");
		exit(1);
	}
	abx = (a+b-1) *(1.-x);
	if ( abx <= 0.8) {
		c2 = (abx * (3. - x) - (b - 1.) *(1. - x))/2.0;
		oup = qgamma(c2,b); //Incomplete Gamma Function Upper Tail
	} else if (abx > 0.8) {
		w1 = pow(b*x, (double) 1./3.);
		w2 = pow(a*(1. - x), (double) 1./3.);
		num = - 3.0 * (w1 * (1. - 1./(9. * b)) - w2 * (1. - 1./ (9. * a)))/sqrt(2);
		den = sqrt((w1 * w1 / b) + (w2 * w2 / a) ); 
		oup = erfc(num/den) / 2.0;
	}
	 
	return oup;
	
}

double ibeta(double x, double a , double b) {
	double oup,temp;
	int ctrlint;
	/*
	 * This program uses two methods described in Abramowitz and Stegun
	 * Handbook.
	 * Method 1. is based on 26.5.8 and uses Continued Fraction. The result
	 * is accepted if the convergence occurs in 1000 iterations.
	 * 
	 * If the program doesn't converge , then the results are discarded and 
	 * method 2. is used.
	 * 
	 * Method 2. uses approximation and is based on equations 26.5.20 (erf)
	 * and 26.5.21 (pgamma/qgamma) 
	 */ 
	if (x < 0. || x > 1.) {
		printf("x only accepts real values between 0.0 and 1.0");
		exit(1);
	}
	if (a < 0. || b < 0.) {
		printf(" The Two Inputs should be nonnegative and real");
		exit(1);
	}
	
	ctrlint = 0;
	
	// Method 1
	
	if ( x < (a+1.)/(a+b+2.)) {
		temp = exp(gamma_log(a+b) - gamma_log(a+1.) - gamma_log(b) + a*log(x) + b*log1p(-x));
		oup = temp * bfrac(x,a,b,&ctrlint); 
	}
	
	if ( x >= (a+1.)/(a+b+2.)) {
		temp = exp(gamma_log(a+b) - gamma_log(b+1.) - gamma_log(a) + a*log(x) + b*log1p(-x));
		oup = 1. - temp * bfrac(1. - x,b,a,&ctrlint); 
	}
	
	// Method 2 is invoked only if ctrlint = 1000
	if (ctrlint == 1000) {
		oup = 0.0;
		oup = ibeta_appx(x,a,b);
	}
	
	return oup;
	
}

double ibetac(double x, double a , double b) {
	double oup;
	
	oup = ibeta(1. - x,b,a);
	
	return oup;
}

double ibetad(double x, double a , double b) {
	double oup;
	/*
	 * This function computes derivative of Incomplete Beta Function
	 *   
	 * 
	 */
	 
	 oup = betapdf(x,a,b);
	  
	return oup;
}

double betapdf(double x, double a , double b) {
	double oup,la,lb;
	/*
	 * The pdf of Beta distribution function is given by
	 * x^(a-1) * (1-x)^(b-1) / B(a,b)
	 * 
	 * or
	 * 
	 * exp((a-1)*log(x) + (b-1)*log(1-x) - beta_log(a,b))
	 * 
	 * Ref 26.1.33 Abramowitz and Stegun
	 */ 
	 
	if (x < 0. || x > 1.) {
		printf("x only accepts real values between 0.0 and 1.0");
		exit(1);
	}
	if (a < 0. || b < 0.) {
		printf(" The Two Inputs should be nonnegative and real");
		exit(1);
	} 
	
	if (a == 1. || x == 0.) {
		la = 0.;
	} else {
		la = (a - 1.) * log(x);
	}
	
	if (b == 1. || x == 1.) {
		lb = 0.;
	} else {
		lb = (b - 1.) * log(1. - x);
	}
	
	oup = exp(la + lb - beta_log(a,b));
	
	return oup;
}

double betacdf(double x, double a , double b) {
	
	double oup;
	oup = ibeta(x,a,b);
	// Since ibeta's method 2 uses approximation it is possible that
	// the cdf function thus calculated may slightly overshoot 1.0 probability
	// value so we check for any discrepancies and correct them.
	
	if (oup > 1.0) {
		oup = 1.0;
	}
	return oup;
	
}

double r8_max ( double x, double y )

/******************************************************************************/
// Unmodified code written by John Burkardt
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 August 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
  }
  return value;
}

double betainv(double alpha, double p, double q) {
	/*
	 *   Author:

    Original FORTRAN77 version by GW Cran, KJ Martin, GE Thomas.
    C version by John Burkardt. 
    * 
    * Slightly modified to fit the rest of the code. (June 17 2013)
     
    Licensing:

    This code is distributed under the GNU LGPL license. 

    Reference:

    GW Cran, KJ Martin, GE Thomas,
    Remark AS R19 and Algorithm AS 109:
    A Remark on Algorithms AS 63: The Incomplete Beta Integral
    and AS 64: Inverse of the Incomplete Beta Integeral,
    Applied Statistics,
    Volume 26, Number 1, 1977, pages 111-114.
	 */ 
	 
  double a;
  double acu;
  double adj;
  double fpu;
  double g;
  double h;
  int iex;
  int indx;
  double pp;
  double prev;
  double qq;
  double r;
  double s;
  double sae = -37.0;
  double sq;
  double t;
  double tx;
  double value;
  double w;
  double xin;
  double y;
  double yprev;
  
  double beta;
 
  beta = beta_log(p,q);

  fpu = pow ( 10.0, sae );
  value = alpha;
  
    if ( p < 0.0 || q < 0.0 )
  {
    printf("parameters should be non-negatve and real");
    exit(1);
    
  }
  
    if ( alpha < 0.0 || 1.0 < alpha )
  {
	printf("probability values must lie between 0 and 1");
    exit(1);
  }
  
    if ( alpha == 0.0 )
  {
    value = 0.0;
    return value;
  }

  if ( alpha == 1.0 )
  {
    value = 1.0;
    return value;
  }
  

  if ( 0.5 < alpha )
  {
    a = 1.0 - alpha;
    pp = q;
    qq = p;
    indx = 1;
  }
  else
  {
    a = alpha;
    pp = p;
    qq = q;
    indx = 0;
  }

// initial approximation

  r = sqrt ( - log ( a * a ) );

  y = r - ( 2.30753 + 0.27061 * r ) 
    / ( 1.0 + ( 0.99229 + 0.04481 * r ) * r );

  if ( 1.0 < pp && 1.0 < qq )
  {
    r = ( y * y - 3.0 ) / 6.0;
    s = 1.0 / ( pp + pp - 1.0 );
    t = 1.0 / ( qq + qq - 1.0 );
    h = 2.0 / ( s + t );
    w = y * sqrt ( h + r ) / h - ( t - s ) 
      * ( r + 5.0 / 6.0 - 2.0 / ( 3.0 * h ) );
    value = pp / ( pp + qq * exp ( w + w ) );
  }
  else
  {
    r = qq + qq;
    t = 1.0 / ( 9.0 * qq );
    t = r * pow ( 1.0 - t + y * sqrt ( t ), 3 );

    if ( t <= 0.0 )
    {
      value = 1.0 - exp ( ( log ( ( 1.0 - a ) * qq ) + beta ) / qq );
    }
    else
    {
      t = ( 4.0 * pp + r - 2.0 ) / t;

      if ( t <= 1.0 )
      {
        value = exp ( ( log ( a * pp ) + beta ) / pp );
      }
      else
      {
        value = 1.0 - 2.0 / ( t + 1.0 );
      }
    }
  }

//Modified Newton-Raphson method
  r = 1.0 - pp;
  t = 1.0 - qq;
  yprev = 0.0;
  sq = 1.0;
  prev = 1.0;

  if ( value < 0.0001 )
  {
    value = 0.0001;
  }

  if ( 0.9999 < value )
  {
    value = 0.9999;
  }

  iex = (int) r8_max ( - 5.0 / pp / pp - 1.0 / pow ( a, 0.2 ) - 13.0, sae );

  acu = pow ( 10.0, iex );

  for ( ; ; )
  {
    y = betacdf( value, pp, qq);

    xin = value;
    y = ( y - a ) * exp ( beta + r * log ( xin ) + t * log ( 1.0 - xin ) );

    if ( y * yprev <= 0.0 )
    {
      prev = r8_max ( sq, fpu );
    }

    g = 1.0;

    for ( ; ; )
    {
      for ( ; ; )
      {
        adj = g * y;
        sq = adj * adj;

        if ( sq < prev )
        {
          tx = value - adj;

          if ( 0.0 <= tx && tx <= 1.0 )
          {
            break;
          }
        }
        g = g / 3.0;
      }

      if ( prev <= acu )
      {
        if ( indx )
        {
          value = 1.0 - value;
        }
        return value;
      }

      if ( y * y <= acu )
      {
        if ( indx )
        {
          value = 1.0 - value;
        }
        return value;
      }

      if ( tx != 0.0 && tx != 1.0 )
      {
        break;
      }

      g = g / 3.0;
    }

    if ( tx == value )
    {
      break;
    }

    value = tx;
    yprev = y;
  }

  if ( indx )
  {
    value = 1.0 - value;
  }

	return value;	
}
