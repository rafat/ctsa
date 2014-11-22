/*
 * nls.c
 *
 *  Created on: May 21, 2014
 *      Author: Rafat Hussain
 */

#include "nls.h"

double enorm(double *x, int N) {
	double enrm;
	int i;
    double agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,x1max,x3max,zero;

	/*
	 * C Translation of Fortran Code
	 * Argonne national laboratory. minpack project. march 1980.
      burton s. garbow, kenneth e. hillstrom, jorge j. more
	 */

    one = 1.0;
    zero = 0.0;
    rdwarf = 3.834e-20;
    rgiant = 1.304e19;
    s1 = s2 = s3 = x1max = x3max = zero;
    floatn = N;
    agiant = rgiant / floatn;

    for(i = 0; i < N;++i) {
    	xabs = fabs(x[i]);
    	if (xabs >= agiant) {
    		//sum for large components
    		if (xabs > x1max) {
    			s1 = one + s1*(x1max/xabs)*(x1max/xabs);
    			x1max = xabs;
    		} else {
    			s1 = s1 + (xabs/x1max)*(xabs/x1max);
    		}
    	} else if (xabs <= rdwarf) {
    		//sum for small components
    		if (xabs > x3max) {
                s3 = one + s3*(x3max/xabs)*(x3max/xabs);
                x3max = xabs;
    		} else {
    			if (xabs != zero) {
    				s3 = s3 + (xabs/x3max)*(xabs/x3max);
    			}
    		}

    	} else {
    		//sum for intermediate components.
    		s2 += xabs*xabs;
    	}
    }

    if (s1 != zero) {
    	enrm = x1max*sqrt(s1+(s2/x1max)/x1max);
    } else {
    	if (s2 != zero) {
    		if (s2 >= x3max) {
    			enrm = sqrt(s2*(one+(x3max/s2)*(x3max*s3)));
    		}
    		if (s2 < x3max) {
    			enrm = sqrt(x3max*((s2/x3max)+(x3max*s3)));
    		}
    	} else {
    		enrm = x3max*sqrt(s3);
    	}
    }

	return enrm;
}

void qrfac(double *A, int M, int N, int lda, int pivot, int *ipvt, int lipvt,double *rdiag, double *acnorm,double eps) {
	int i,j,jp1,k,kmax,minmn,t;
    double ajnorm,epsmch,one,p05,sum,temp,zero,temp2,pmaxval;
    double *AT,*wa,*wa2;

 /*
     * This routine is a C translation of Fortran Code by
     *     argonne national laboratory. minpack project. march 1980.
      burton s. garbow, kenneth e. hillstrom, jorge j. more
     *
     * M is a positive integer input variable set to the number
         of rows of a.

       N is a positive integer input variable set to the number
         of columns of a.

       A is an M by N array. on input a contains the matrix for
         which the qr factorization is to be computed. on output
         the strict upper trapezoidal part of a contains the strict
         upper trapezoidal part of r, and the lower trapezoidal
         part of a contains a factored form of q (the non-trivial
         elements of the u vectors described above).

       lda is a positive integer input variable not less than m
         which specifies the leading dimension of the array a.

       pivot is an integer input variable. if pivot is set to 1,
         then column pivoting is enforced. if pivot is set to 0,
         then no column pivoting is done.

       ipvt is an integer output array of length lipvt. ipvt
         defines the permutation matrix p such that a*p = q*r.
         column j of p is column ipvt(j) of the identity matrix.
         if pivot is false, ipvt is not referenced.

       lipvt is a positive integer input variable. if pivot is false,
         then lipvt may be as small as 1. if pivot is true, then
         lipvt must be at least n.

       rdiag is an output array of length N which contains the
         diagonal elements of r.

       acnorm is an output array of length N which contains the
         norms of the corresponding columns of the input matrix a.
         if this information is not needed, then acnorm can coincide
         with rdiag.
     *
     */


    if (pivot != 0 && pivot != 1) {
    	printf("Pivot only takes binary values 0 and 1 \n");
    	exit(-1);
    }

    AT = (double*) malloc(sizeof(double) *N*M);
    wa = (double*) malloc(sizeof(double) *N);
    wa2 = (double*) malloc(sizeof(double) *M);

    one = 1.0; zero = 0.0; p05 = 5.0e-02;
    epsmch = eps;

    mtranspose(A,M,N,AT);// AT is size NXM

    //compute the initial column norms and initialize several arrays.

    for(j = 0; j < N;++j) {
    	acnorm[j] = enorm(AT+j*M,M);
    	rdiag[j] = acnorm[j];
    	wa[j] = rdiag[j];
    	if (pivot == 1) {
    		ipvt[j] = j;
    	}
    }

    //reduce a to r with householder transformations.

    if (M < N) {
    	minmn = M;
    } else {
    	minmn = N;
    }

    for (j = 0; j < minmn;++j) {
    	if (pivot == 1) {
    		//bring the column of largest norm into the pivot position.
    		kmax = j;
    		for(k = j; k < N;++k) {
    			if (rdiag[k] > rdiag[kmax]) {
    				kmax = k;
    			}
    		}
			if (kmax != j) {
				for(i = 0; i < M;++i) {
					t = i * N;
					temp = A[t+j];
					A[t+j] = A[t+kmax];
					A[t+kmax] = temp;
				}
				rdiag[kmax] = rdiag[j];
				wa[kmax] = wa[j];
				k = ipvt[j];
				ipvt[j] = ipvt[kmax];
				ipvt[kmax] = k;
			}
    	}
    	//        compute the householder transformation to reduce the
    	//       j-th column of a to a multiple of the j-th unit vector.
    	t = j * N + j;

    	for(i = 0; i < M-j;++i) {
    		wa2[i] = A[t+i*N];
    	}
    	ajnorm = enorm(wa2,M-j);
    	if (ajnorm != zero) {
    		if (A[t] < zero) {
    			ajnorm = - ajnorm;
    		}
    		for(i = j; i < M;++i) {
    			A[i*N+j] /= ajnorm;
    		}
    		A[t] += one;
    		//        apply the transformation to the remaining columns
    		//        and update the norms.

    		jp1 = j + 1; // Breakpoint
    		if (N >= jp1+1) {
    			for(k = jp1; k < N;++k) {
    				sum = zero;
    				for(i = j; i < M;++i) {
    					sum += (A[i*N+j] * A[i*N+k]);
    				}
    				temp = sum / A[t];
    				for(i = j; i < M;++i) {
    					A[i*N+k] -= (temp * A[i*N+j]);
    				}
    				// Breakpoint 2
    				if (pivot == 1 && rdiag[k] != zero) {
    					temp = A[j*N+k] / rdiag[k];
						pmaxval = pmax(zero, one - temp*temp);
    					rdiag[k] = rdiag[k]*sqrt(pmaxval);
    					temp2 = (p05*(rdiag[k]/wa[k]));
    					temp2 = temp2 * temp2;
    					if (temp2 <= epsmch) {
    				    	for(i = 0; i < M-j-1;++i) {
    				    		wa2[i] = A[jp1*N+k+i*N];
    				    	}
    				    	rdiag[k] = enorm(wa2,M-j-1);
    				    	wa[k] = rdiag[k];
    					}
    				}
    			}
    		}
    	}
    	rdiag[j] = -ajnorm;
    }

    free(AT);
    free(wa);
    free(wa2);
}

void qrsolv(double *r,int ldr,int N,int *ipvt,double *diag,double *qtb,double *x,double *sdiag) {
	int i,j,jp1,k,kp1,l,nsing,t;
	double cos,cotan,p5,p25,qtbpj,sin,sum,tan,temp,zero;
	double *wa;

	/*
	 *   This routine is a C translation of Fortran Code by
     *     argonne national laboratory. minpack project. march 1980.
      	  burton s. garbow, kenneth e. hillstrom, jorge j. more
	 * N is a positive integer input variable set to the order of r.

       r is an N by N array. on input the full upper triangle
         must contain the full upper triangle of the matrix r.
         on output the full upper triangle is unaltered, and the
         strict lower triangle contains the strict upper triangle
         (transposed) of the upper triangular matrix s.

       ldr is a positive integer input variable not less than n
         which specifies the leading dimension of the array r.

       ipvt is an integer input array of length N which defines the
         permutation matrix p such that a*p = q*r. column j of p
         is column ipvt(j) of the identity matrix.

       diag is an input array of length N which must contain the
         diagonal elements of the matrix d.

       qtb is an input array of length N which must contain the first
         N elements of the vector (q transpose)*b.

       x is an output array of length N which contains the least
         squares solution of the system a*x = b, d*x = 0.

       sdiag is an output array of length N which contains the
         diagonal elements of the upper triangular matrix s.
	 */

	wa = (double*) malloc(sizeof(double) *N);

	p5 = 5.0e-1;
	p25 = 2.5e-1;
	zero = 0.0;

	//     copy r and (q transpose)*b to preserve input and initialize s.
	//     in particular, save the diagonal elements of r in x.

	for(j = 0; j < N;++j) {
		for(i = j; i < N;++i) {
			r[i*N+j] = r[j*N+i];
		}
		x[j] = r[j*N+j];
		wa[j] = qtb[j];
	}

	//eliminate the diagonal matrix d using a givens rotation.

	for(j = 0; j < N;++j) { //100

		//        prepare the row of d to be eliminated, locating the
		//        diagonal element using p from the qr factorization.

		l = ipvt[j];

		if (diag[l] != zero) { //90
			for(k = j; k < N;++k) {
				sdiag[k] = zero;
			}

			sdiag[j] = diag[l];

			//        the transformations to eliminate the row of d
			//        modify only a single element of (q transpose)*b
			//        beyond the first n, which is initially zero.

			qtbpj = zero;

			for(k = j; k < N;++k) { //80
				t = k * N;
				//     determine a givens rotation which eliminates the
				//     appropriate element in the current row of d.

				if (sdiag[k] != zero) { //70
					if (fabs(r[t+k]) < fabs(sdiag[k])) {
						cotan = r[t+k]/sdiag[k];
			            sin = p5/sqrt(p25+p25*cotan*cotan);
			            cos = sin*cotan;
					} else {
			            tan = sdiag[k]/r[t+k];
			            cos = p5/sqrt(p25+p25*tan*tan);
			            sin = cos*tan;
					}
					//           compute the modified diagonal element of r and
					//           the modified element of ((q transpose)*b,0).
		            r[t+k] = cos*r[t+k] + sin*sdiag[k];
		            temp = cos*wa[k] + sin*qtbpj;
		            qtbpj = -sin*wa[k] + cos*qtbpj;
		            wa[k] = temp;

		            //           accumulate the tranformation in the row of s.
		            kp1 = k + 1;

		            if (N >= kp1+1) {//71
		            	for(i = kp1;i < N;++i) {
		                    temp = cos*r[i*N+k] + sin*sdiag[i];
		                    sdiag[i] = -sin*r[i*N+k] + cos*sdiag[i];
		                    r[i*N+k] = temp;
		            	}
		            }//71
				}//70
			}//80
		} //90
		//        store the diagonal element of s and restore
		//        the corresponding diagonal element of r.

        sdiag[j] = r[j*N+j];
        r[j*N+j] = x[j];

	}//100

	//     solve the triangular system for z. if the system is
	//     singular, then obtain a least squares solution.

	nsing = N;

	for(j = 1; j <= N;++j) {
		if (sdiag[j-1] == zero && nsing == N) {
			nsing = j - 1;
		}
		if (nsing < N) {
			wa[j-1] = zero;
		}
	}

	if (nsing >= 1) {//150
		for(k = 1; k <= nsing;++k) {
	         j = nsing - k + 1;
	         sum = zero;
	         jp1 = j + 1;

	         if (nsing >= jp1) {
	        	 for (i = jp1; i <= nsing;++i) {
	        		 sum = sum + r[(i-1)*N+j-1]*wa[i-1];
	        	 }
	         }
	         wa[j-1] = (wa[j-1] - sum)/sdiag[j-1];
		}
	}//150

	//     permute the components of z back to components of x.

	for (j = 0; j < N;++j) {
        l = ipvt[j];
        x[l] = wa[j];
	}


	free(wa);
}

void fdjac2(custom_funcmult *funcmult, double *x, int M, int N, double *fvec, double *fjac, int ldfjac,
		double epsfcn,double eps) {
	int i,j;
	double epsmch,h,temp,zero;
	double *wa;

	zero = 0.0;
	epsmch = eps;
	eps = sqrt(pmax(epsfcn,epsmch));

	wa = (double*) malloc(sizeof(double) *M);

	for(j = 0; j < N;++j) {
        temp = x[j];
        h = eps*fabs(temp);
        if (h == zero) {
        	h = eps;
        }
        x[j] = temp + h;
        FUNCMULT_EVAL(funcmult,x,M,N,wa);
        x[j] = temp;
        for(i = 0; i < M;++i) {
        	fjac[i*N+j] = (wa[i] - fvec[i])/h;
        }
	}

	free(wa);

}

void lmpar(double *r,int ldr,int N,int *ipvt,double *diag,double *qtb,double delta,double *par,double *x,double *sdiag) {
	int i,iter,j,jm1,jp1,k,l,nsing;
	double dxnorm,dwarf,fp,gnorm,parc,parl,paru,p1,p001,sum,temp,zero;
	double *wa1,*wa2;
	/*
	*   This routine is a C translation of Fortran Code by
    *     argonne national laboratory. minpack project. march 1980.
     	  burton s. garbow, kenneth e. hillstrom, jorge j. more

	 * N is a positive integer input variable set to the order of r.

       r is an N by N array. on input the full upper triangle
         must contain the full upper triangle of the matrix r.
         on output the full upper triangle is unaltered, and the
         strict lower triangle contains the strict upper triangle
         (transposed) of the upper triangular matrix s.

       ldr is a positive integer input variable not less than n
         which specifies the leading dimension of the array r.

       ipvt is an integer input array of length N which defines the
         permutation matrix p such that a*p = q*r. column j of p
         is column ipvt(j) of the identity matrix.

       diag is an input array of length N which must contain the
         diagonal elements of the matrix d.

       qtb is an input array of length N which must contain the first
         N elements of the vector (q transpose)*b.

       delta is a positive input variable which specifies an upper
         bound on the euclidean norm of d*x.

       par is a nonnegative variable. on input par contains an
         initial estimate of the levenberg-marquardt parameter.
         on output par contains the final estimate.

       x is an output array of length N which contains the least
         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
         for the output par.

       sdiag is an output array of length N which contains the
         diagonal elements of the upper triangular matrix s.
    */

	wa1 = (double*) malloc(sizeof(double) *N);
	wa2 = (double*) malloc(sizeof(double) *N);

	p1 = 1.0e-01;
	p001 = 1.0e-03;
	zero = 0.0;
	dwarf = 2.22507385852e-308;

	//     compute and store in x the gauss-newton direction. if the
	//     jacobian is rank-deficient, obtain a least squares solution.

	nsing = N;

	for(j = 1; j <= N;++j) {
		wa1[j-1] = qtb[j-1];
		if (r[(j-1)*N+j-1] == zero && nsing == N) {
			nsing = j - 1;
		}
		if (nsing < N) {
			wa1[j-1] = zero;
		}
	}

	if (nsing >= 1) {//50
		for(k = 1; k <= nsing;++k) {
	         j = nsing - k + 1;
	         wa1[j-1] = wa1[j-1]/r[(j-1)*N+j-1];
	         temp = wa1[j-1];
	         jm1 = j - 1;
	         if (jm1 >= 1) {
	        	 for(i = 1; i <= jm1;++i) {
	        		 wa1[i-1] = wa1[i-1] - r[(i-1)*N+j-1]*temp;
	        	 }
	         }
		}
	}//50

	for (j = 0; j < N;++j) {
        l = ipvt[j];
        x[l] = wa1[j];
	}

	//     initialize the iteration counter.
	//     evaluate the function at the origin, and test
	//     for acceptance of the gauss-newton direction.

	iter = 0;

	for(j = 0; j < N;++j) {
		wa2[j] = diag[j]*x[j];
	}

	dxnorm = enorm(wa2,N);
	fp = dxnorm - delta;

	if (fp > p1*delta) {//220
		//     if the jacobian is not rank deficient, the newton
		//    step provides a lower bound, parl, for the zero of
		//     the function. otherwise set this bound to zero.
		parl = zero;

		if (nsing >= N) { //120 nsing only takes values upto N
			for(j = 0; j < N;++j) {
				l = ipvt[j];
				wa1[j] = diag[l]*(wa2[l]/dxnorm);
			}

			for(j = 0; j < N;++j) {//110
		         sum = zero;
		         jm1 = j - 1;
		         if (jm1 >= 0) {//100
		        	 for(i = 0; i <= jm1;++i) {//check
		        		 sum = sum + r[i*N+j]*wa1[i];
		        	 }
		         }//100
		         wa1[j] = (wa1[j] - sum)/r[j*N+j];
			}//110
		      temp = enorm(wa1,N);
		      parl = ((fp/delta)/temp)/temp;
		}//120

		//     calculate an upper bound, paru, for the zero of the function.

		for(j = 0; j < N;++j) {//140
			sum = zero;
			for(i = 0; i <= j;++i) {//check
				sum = sum + r[i*N+j]*qtb[i];
			}
	         l = ipvt[j];
	         wa1[j] = sum/diag[l];
		}//140
		gnorm = enorm(wa1,N);
		paru = gnorm/delta;

		if (paru == zero) {
			paru = dwarf/pmin(delta,p1);
		}

		//     if the input par lies outside of the interval (parl,paru),
		//     set par to the closer endpoint.

	      *par = pmax(*par,parl);
	      *par = pmin(*par,paru);

	      if (*par == zero) {
	    	  *par = gnorm/dxnorm;
	      }
	      //Iteration begins
	      while(1) {
	    	  iter++;
	    	  //        evaluate the function at the current value of par.
	    	  if (*par == zero) {
	    		  *par = pmax(dwarf,p001*paru);
	    	  }
	    	  temp = sqrt(*par);
	    	  for(j = 0; j < N;++j) {
	    		  wa1[j] = temp*diag[j];
	    	  }

	    	  qrsolv(r,ldr,N,ipvt,wa1,qtb,x,sdiag);
	    	  for(j = 0; j < N;++j) {
	              wa2[j] = diag[j]*x[j];
	    	  }

	          dxnorm = enorm(wa2,N);
	          temp = fp;
	          fp = dxnorm - delta;

	          //        if the function is small enough, accept the current value
	          //        of par. also test for the exceptional cases where parl
	          //        is zero or the number of iterations has reached 10.
	          if (fabs(fp) <= p1*delta) {
	        	  break;
	          }

	          if (iter == 10) {
	        	  break;
	          }

	          if (parl == zero && fp <= temp && temp < zero) {
	        	  break;
	          }

	          //        compute the newton correction.

	          for(j = 0;j < N;++j) {//180
	              l = ipvt[j];
	              wa1[j] = diag[l]*(wa2[l]/dxnorm);
	          }//180

	          for(j = 0; j < N;++j) {//210
	              wa1[j] = wa1[j]/sdiag[j];
	              temp = wa1[j];
	              jp1 = j + 1;
	              if (N >= jp1+1) {
	            	  for(i = jp1; i < N;++i) {
	            		  wa1[i] = wa1[i] - r[i*N+j]*temp;
	            	  }
	              }
	          }//210
	          temp = enorm(wa1,N);
	          parc = ((fp/delta)/temp)/temp;
	          //        depending on the sign of the function, update parl or paru.

	          if (fp > zero) {
	        	  parl = pmax(parl,*par);
	          }
	          if (fp < zero) {
	        	  paru = pmin(paru,*par);
	          }

	          //       compute an improved estimate for par.
	          *par = pmax(parl,*par+parc);

	      }

	}//220

	if (iter == 0) {
		*par = zero;
	}

	free(wa1);
	free(wa2);
}

int lmder(custom_funcmult *funcmult, custom_jacobian *jacobian, double *x, int M, int N,
		double *fvec,double *fjac,int ldfjac,int maxfev,double *diag,int mode,double factor,int nprint,
		double eps,double ftol,double gtol,double xtol,int *nfev,int *njev,int *ipvt, double *qtf) {
	int info;
	int i,j,l,iter;
    double actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,
    sum,temp,temp1,temp2,xnorm,zero;
    double *wa1,*wa2,*wa3,*wa4;

    /*
     * 	*   This routine is a C translation of Fortran Code by
    *     argonne national laboratory. minpack project. march 1980.
     	  burton s. garbow, kenneth e. hillstrom, jorge j. more
     *  M is a positive integer input variable set to the number
c         of functions.
c
c       N is a positive integer input variable set to the number
c         of variables. N must not exceed M.
c
c       x is an array of length N. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length M which contains
c         the functions evaluated at the output x.
c
c       fjac is an output M by N array. the upper N by N submatrix
c         of fjac contains an upper triangular matrix r with
c         diagonal elements of nonincreasing magnitude such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final
c         calculated jacobian. column j of p is column ipvt(j)
c         (see below) of the identity matrix. the lower trapezoidal
c         part of fjac contains information generated during
c         the computation of r.
c
c       ldfjac is a positive integer input variable not less than M
c         which specifies the leading dimension of the array fjac.
c
c       ftol is a nonnegative input variable. termination
c         occurs when both the actual and predicted relative
c         reductions in the sum of squares are at most ftol.
c         therefore, ftol measures the relative error desired
c         in the sum of squares.
c
c       xtol is a nonnegative input variable. termination
c         occurs when the relative error between two consecutive
c         iterates is at most xtol. therefore, xtol measures the
c         relative error desired in the approximate solution.
c
c       gtol is a nonnegative input variable. termination
c         occurs when the cosine of the angle between fvec and
c         any column of the jacobian is at most gtol in absolute
c         value. therefore, gtol measures the orthogonality
c         desired between the function vector and the columns
c         of the jacobian.
c
c       maxfev is a positive integer input variable. termination
c         occurs when the number of calls to fcn with iflag = 1
c         has reached maxfev.
c
c       diag is an array of length N. if mode = 1 (see
c         below), diag is internally set. if mode = 2, diag
c         must contain positive entries that serve as
c         multiplicative scale factors for the variables.
c
c       mode is an integer input variable. if mode = 1, the
c         variables will be scaled internally. if mode = 2,
c         the scaling is specified by the input diag. other
c         values of mode are equivalent to mode = 1.
c
c       factor is a positive input variable used in determining the
c         initial step bound. this bound is set to the product of
c         factor and the euclidean norm of diag*x if nonzero, or else
c         to factor itself. in most cases factor should lie in the
c         interval (.1,100.).100. is a generally recommended value.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive. in this case,
c         fcn is called with iflag = 0 at the beginning of the first
c         iteration and every nprint iterations thereafter and
c         immediately prior to return, with x, fvec, and fjac
c         available for printing. fvec and fjac should not be
c         altered. if nprint is not positive, no special calls
c         of fcn with iflag = 0 are made.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  both actual and predicted relative reductions
c                   in the sum of squares are at most ftol.
c
c         info = 2  relative error between two consecutive iterates
c                   is at most xtol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  the cosine of the angle between fvec and any
c                   column of the jacobian is at most gtol in
c                   absolute value.
c
c         info = 5  number of calls to fcn with iflag = 1 has
c                   reached maxfev.
c
c         info = 6  ftol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  xtol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c         info = 8  gtol is too small. fvec is orthogonal to the
c                   columns of the jacobian to machine precision.
c
c       nfev is an integer output variable set to the number of
c         calls to fcn with iflag = 1.
c
c       njev is an integer output variable set to the number of
c         calls to fcn with iflag = 2.
c
c       ipvt is an integer output array of length N. ipvt
c         defines a permutation matrix p such that jac*p = q*r,
c         where jac is the final calculated jacobian, q is
c         orthogonal (not stored), and r is upper triangular
c         with diagonal elements of nonincreasing magnitude.
c         column j of p is column ipvt(j) of the identity matrix.
c
c       qtf is an output array of length N which contains
c         the first n elements of the vector (q transpose)*fvec.
     */

	wa1 = (double*) malloc(sizeof(double) *N);
	wa2 = (double*) malloc(sizeof(double) *N);
	wa3 = (double*) malloc(sizeof(double) *N);
	wa4 = (double*) malloc(sizeof(double) *M);

    one = 1.0;
    zero = 0.0;
    p1 = 1.0e-1; p5 = 5.0e-1; p25 = 2.5e-1; p75 = 7.5e-1; p0001 = 1.0e-4;
    epsmch = eps;

    info = 0;
    *nfev = 0;
    *njev = 0;

    if (N <= 0 || M < N || ldfjac < M || ftol < zero || xtol < zero || gtol < zero || maxfev <= 0 || factor <= zero) {
    	return info;
    }
    if (mode == 2) {
		for(j = 0; j < N;++j) {
			if (diag[j] <= 0.0) {
				return info;
			}
		}
    }

    //     evaluate the function at the starting point
    //     and calculate its norm.

    FUNCMULT_EVAL(funcmult,x,M,N,fvec);
    *nfev= 1;
    fnorm = enorm(fvec,M);

    //     initialize levenberg-marquardt parameter and iteration counter.
    par = zero;
    iter = 1;
    ratio = zero;

    //     beginning of the outer loop.

    while(1) {
    	//        calculate the jacobian matrix.
    	ratio = zero;
    	JACOBIAN_EVAL(jacobian,x,M,N,fjac);
    	*njev = *njev +1;

    	//        compute the qr factorization of the jacobian.

    	qrfac(fjac,M,N,ldfjac,1,ipvt,N,wa1,wa2,eps);

    	//        on the first iteration and if mode is 1, scale according
    	//        to the norms of the columns of the initial jacobian.

    	if (iter == 1) {//80
    		if (mode != 2) {//60
    			for(j = 0;j < N;++j) {
    				diag[j] = wa2[j];
    				if (wa2[j] == zero) {
    					diag[j] = one;
    				}
    			}
    		}//60

    		//        on the first iteration, calculate the norm of the scaled x
    		//        and initialize the step bound delta.

    		for(j = 0; j < N;++j) {
    			wa3[j] = diag[j]*x[j];
    		}
            xnorm = enorm(wa3,N);
            delta = factor*xnorm;

            if (delta == zero) {
            	delta = factor;
            }

    	}//80

        //        form (q transpose)*fvec and store the first n components in
        //        qtf.

    	for(i = 0; i < M;++i) {
    		wa4[i] = fvec[i];
    	}

    	for(j = 0; j < N;++j) {//130
    		if (fjac[j*N+j] != zero) {//120
    			sum = zero;
    			for(i = j; i < M;++i) {//100
    				sum = sum + fjac[i*N+j]*wa4[i];
    			}//100
    			temp = -sum/fjac[j*N+j];
    			for(i = j; i < M;++i) {//110
    				wa4[i] = wa4[i] + fjac[i*N+j]*temp;
    			}//110
    		}//120
            fjac[j*N+j] = wa1[j];
            qtf[j] = wa4[j];
    	}//130

    	//        compute the norm of the scaled gradient.
    	gnorm = zero;

    	if (fnorm != zero) {//170
    		for(j = 0; j < N;++j) {//160
    			l = ipvt[j];
    			if (wa2[l] != zero) {//150
    				sum = zero;
    				for(i = 0; i <= j;++i) { //140
    					sum = sum + fjac[i*N+j]*(qtf[i]/fnorm);
    				}//140
    				gnorm = pmax(gnorm,fabs(sum/wa2[l]));
    			}//150
    		}//160
    	}//170

    	//        test for convergence of the gradient norm.
    	if (gnorm <= gtol) {
    		info = 4;
    	}
    	if (info != 0) {
    		break;
    	}

    	//        rescale if necessary.
    	if (mode != 2) { //190
    		for(j = 0; j < N;++j) {
    			diag[j] = pmax(diag[j],wa2[j]);
    		}
    	}//190

    	//        beginning of the inner loop.

    	while(ratio < p0001) {
    		//           determine the levenberg-marquardt parameter.
    		lmpar(fjac,ldfjac,N,ipvt,diag,qtf,delta,&par,wa1,wa2);
    		//           store the direction p and x + p. calculate the norm of p.
    		for(j = 0; j < N;++j) {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j]*wa1[j];
    		}
    		pnorm = enorm(wa3,N);
    		//           on the first iteration, adjust the initial step bound.
    		if (iter == 1) {
    			delta = pmin(delta,pnorm);
    		}
    		//           evaluate the function at x + p and calculate its norm.

    		FUNCMULT_EVAL(funcmult,wa2,M,N,wa4);
    		*nfev = *nfev + 1;
    		fnorm1 = enorm(wa4,M);

    		//           compute the scaled actual reduction.

            actred = -one;
            if (p1*fnorm1 < fnorm) {
            	actred = one - (fnorm1/fnorm)*(fnorm1/fnorm);
            }

            //           compute the scaled predicted reduction and
            //           the scaled directional derivative.

            for(j = 0; j < N;++j) {
                wa3[j] = zero;
                l = ipvt[j];
                temp = wa1[l];
                for(i = 0;i <= j;++i) {
                	wa3[i] = wa3[i] + fjac[i*N+j]*temp;
                }
            }

            temp1 = enorm(wa3,N);
            temp1 = temp1/fnorm;
            temp2 = (sqrt(par)*pnorm)/fnorm;
            prered = temp1*temp1 + temp2*temp2/p5;
            dirder = -(temp1*temp1 + temp2*temp2);
            //           compute the ratio of the actual to the predicted
            //           reduction.
            ratio = zero;
            if (prered != zero) {
            	ratio = actred/prered;
            }
            //           update the step bound.

            if (ratio <= p25) {//240
            	if (actred >= zero) {
            		temp = p5;
            	}
            	if (actred < zero) {
            		temp = p5*dirder/(dirder + p5*actred);
            	}
            	if (p1*fnorm1 >= fnorm || temp < p1) {
            		temp = p1;
            	}
                delta = temp*pmin(delta,pnorm/p1);
                par = par/temp;
            } else if (par == zero || ratio >= p75){//240 - 260
                delta = pnorm/p5;
                par = p5*par;
            }//260

            //           test for successful iteration.

            if (ratio >= p0001) {//290
            	//           successful iteration. update x, fvec, and their norms.
            	for(j = 0; j < N;++j) {
                    x[j] = wa2[j];
                    wa2[j] = diag[j]*x[j];
            	}
            	for(i = 0; i < M;++i) {
            		fvec[i] = wa4[i];
            	}
                xnorm = enorm(wa2,N);
                fnorm = fnorm1;
                iter = iter + 1;
            }//290
            //           tests for convergence.
            if ((fabs(actred) <= ftol) && (prered <= ftol) && (p5*ratio <= one)) {
            	info = 1;
            }
            if (delta <= xtol*xnorm) {
            	info = 2;
            }
            if ((fabs(actred) <= ftol) && (prered <= ftol) && (p5*ratio <= one) && (info == 2)) {
            	info = 3;
            }
            if (info != 0) {
            	break;
            }

            //           tests for termination and stringent tolerances.
            if (*nfev >= maxfev) {
            	info = 5;
            }
            if ((fabs(actred) <= epsmch) && (prered <= epsmch) && (p5*ratio <= one)) {
            	info = 6;
            }
            if (delta <= epsmch*xnorm) {
            	info = 7;
            }
            if (gnorm <= epsmch) {
            	info = 8;
            }
            if (info != 0) {
            	break;
            }

    	}

        if (info != 0) {
        	break;
        }


    }


    free(wa1);
    free(wa2);
    free(wa3);
    free(wa4);

	return info;
}

int lmdif(custom_funcmult *funcmult, double *x, int M, int N, double *fvec, double *fjac, int ldfjac,
		int maxfev,double *diag,int mode,double factor,int nprint,double eps,double epsfcn,double ftol,double gtol,
		double xtol,int *nfev,int *njev,int *ipvt, double *qtf) {
	int info;
	int i,j,l,iter;
    double actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,
    sum,temp,temp1,temp2,xnorm,zero;
    double *wa1,*wa2,*wa3,*wa4;

    /*
     * 	*   This routine is a C translation of Fortran Code by
    *     argonne national laboratory. minpack project. march 1980.
     	  burton s. garbow, kenneth e. hillstrom, jorge j. more
     *  M is a positive integer input variable set to the number
c         of functions.
c
c       N is a positive integer input variable set to the number
c         of variables. N must not exceed M.
c
c       x is an array of length N. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length M which contains
c         the functions evaluated at the output x.
c
c       fjac is an output M by N array. the upper N by N submatrix
c         of fjac contains an upper triangular matrix r with
c         diagonal elements of nonincreasing magnitude such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final
c         calculated jacobian. column j of p is column ipvt(j)
c         (see below) of the identity matrix. the lower trapezoidal
c         part of fjac contains information generated during
c         the computation of r.
c
c       ldfjac is a positive integer input variable not less than M
c         which specifies the leading dimension of the array fjac.
c
c       ftol is a nonnegative input variable. termination
c         occurs when both the actual and predicted relative
c         reductions in the sum of squares are at most ftol.
c         therefore, ftol measures the relative error desired
c         in the sum of squares.
c
c       xtol is a nonnegative input variable. termination
c         occurs when the relative error between two consecutive
c         iterates is at most xtol. therefore, xtol measures the
c         relative error desired in the approximate solution.
c
c       gtol is a nonnegative input variable. termination
c         occurs when the cosine of the angle between fvec and
c         any column of the jacobian is at most gtol in absolute
c         value. therefore, gtol measures the orthogonality
c         desired between the function vector and the columns
c         of the jacobian.
c
c       maxfev is a positive integer input variable. termination
c         occurs when the number of calls to fcn with iflag = 1
c         has reached maxfev.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. if epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       diag is an array of length N. if mode = 1 (see
c         below), diag is internally set. if mode = 2, diag
c         must contain positive entries that serve as
c         multiplicative scale factors for the variables.
c
c       mode is an integer input variable. if mode = 1, the
c         variables will be scaled internally. if mode = 2,
c         the scaling is specified by the input diag. other
c         values of mode are equivalent to mode = 1.
c
c       factor is a positive input variable used in determining the
c         initial step bound. this bound is set to the product of
c         factor and the euclidean norm of diag*x if nonzero, or else
c         to factor itself. in most cases factor should lie in the
c         interval (.1,100.).100. is a generally recommended value.
c
c       nprint is an integer input variable that enables controlled
c         printing of iterates if it is positive. in this case,
c         fcn is called with iflag = 0 at the beginning of the first
c         iteration and every nprint iterations thereafter and
c         immediately prior to return, with x, fvec, and fjac
c         available for printing. fvec and fjac should not be
c         altered. if nprint is not positive, no special calls
c         of fcn with iflag = 0 are made.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  both actual and predicted relative reductions
c                   in the sum of squares are at most ftol.
c
c         info = 2  relative error between two consecutive iterates
c                   is at most xtol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  the cosine of the angle between fvec and any
c                   column of the jacobian is at most gtol in
c                   absolute value.
c
c         info = 5  number of calls to fcn with iflag = 1 has
c                   reached maxfev.
c
c         info = 6  ftol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  xtol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c         info = 8  gtol is too small. fvec is orthogonal to the
c                   columns of the jacobian to machine precision.
c
c       nfev is an integer output variable set to the number of
c         calls to fcn with iflag = 1.
c
c       njev is an integer output variable set to the number of
c         calls to fcn with iflag = 2.
c
c       ipvt is an integer output array of length N. ipvt
c         defines a permutation matrix p such that jac*p = q*r,
c         where jac is the final calculated jacobian, q is
c         orthogonal (not stored), and r is upper triangular
c         with diagonal elements of nonincreasing magnitude.
c         column j of p is column ipvt(j) of the identity matrix.
c
c       qtf is an output array of length N which contains
c         the first n elements of the vector (q transpose)*fvec.
     */

	wa1 = (double*) malloc(sizeof(double) *N);
	wa2 = (double*) malloc(sizeof(double) *N);
	wa3 = (double*) malloc(sizeof(double) *N);
	wa4 = (double*) malloc(sizeof(double) *M);

    one = 1.0;
    zero = 0.0;
    p1 = 1.0e-1; p5 = 5.0e-1; p25 = 2.5e-1; p75 = 7.5e-1; p0001 = 1.0e-4;
    epsmch = eps;

    info = 0;
    *nfev = 0;
    *njev = 0;

    if (N <= 0 || M < N || ldfjac < M || ftol < zero || xtol < zero || gtol < zero || maxfev <= 0 || factor <= zero) {
    	return info;
    }

    if (mode == 2) {
		for(j = 0; j < N;++j) {
			if (diag[j] <= 0.0) {
				return info;
			}
		}
    }

    //     evaluate the function at the starting point
    //     and calculate its norm.

    FUNCMULT_EVAL(funcmult,x,M,N,fvec);
    *nfev= 1;
    fnorm = enorm(fvec,M);

    //     initialize levenberg-marquardt parameter and iteration counter.
    par = zero;
    iter = 1;
    ratio = zero;

    //     beginning of the outer loop.

    while(1) {
    	//        calculate the jacobian matrix.
    	ratio = zero;
    	fdjac2(funcmult,x,M,N,fvec,fjac,ldfjac,epsfcn,epsmch);
    	*njev = *njev + N;

    	//        compute the qr factorization of the jacobian.

    	qrfac(fjac,M,N,ldfjac,1,ipvt,N,wa1,wa2,eps);

    	//        on the first iteration and if mode is 1, scale according
    	//        to the norms of the columns of the initial jacobian.

    	if (iter == 1) {//80
    		if (mode != 2) {//60
    			for(j = 0;j < N;++j) {
    				diag[j] = wa2[j];
    				if (wa2[j] == zero) {
    					diag[j] = one;
    				}
    			}
    		}//60

    		//        on the first iteration, calculate the norm of the scaled x
    		//        and initialize the step bound delta.

    		for(j = 0; j < N;++j) {
    			wa3[j] = diag[j]*x[j];
    		}
            xnorm = enorm(wa3,N);
            delta = factor*xnorm;

            if (delta == zero) {
            	delta = factor;
            }

    	}//80

        //        form (q transpose)*fvec and store the first n components in
        //        qtf.

    	for(i = 0; i < M;++i) {
    		wa4[i] = fvec[i];
    	}

    	for(j = 0; j < N;++j) {//130
    		if (fjac[j*N+j] != zero) {//120
    			sum = zero;
    			for(i = j; i < M;++i) {//100
    				sum = sum + fjac[i*N+j]*wa4[i];
    			}//100
    			temp = -sum/fjac[j*N+j];
    			for(i = j; i < M;++i) {//110
    				wa4[i] = wa4[i] + fjac[i*N+j]*temp;
    			}//110
    		}//120
            fjac[j*N+j] = wa1[j];
            qtf[j] = wa4[j];
    	}//130

    	//        compute the norm of the scaled gradient.
    	gnorm = zero;

    	if (fnorm != zero) {//170
    		for(j = 0; j < N;++j) {//160
    			l = ipvt[j];
    			if (wa2[l] != zero) {//150
    				sum = zero;
    				for(i = 0; i <= j;++i) { //140
    					sum = sum + fjac[i*N+j]*(qtf[i]/fnorm);
    				}//140
    				gnorm = pmax(gnorm,fabs(sum/wa2[l]));
    			}//150
    		}//160
    	}//170

    	//        test for convergence of the gradient norm.
    	if (gnorm <= gtol) {
    		info = 4;
    	}
    	if (info != 0) {
    		break;
    	}

    	//        rescale if necessary.
    	if (mode != 2) { //190
    		for(j = 0; j < N;++j) {
    			diag[j] = pmax(diag[j],wa2[j]);
    		}
    	}//190

    	//        beginning of the inner loop.

    	while(ratio < p0001) {
    		//           determine the levenberg-marquardt parameter.
    		lmpar(fjac,ldfjac,N,ipvt,diag,qtf,delta,&par,wa1,wa2);
    		//           store the direction p and x + p. calculate the norm of p.
    		for(j = 0; j < N;++j) {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j]*wa1[j];
    		}
			
    		pnorm = enorm(wa3,N);
    		//           on the first iteration, adjust the initial step bound.
    		if (iter == 1) {
    			delta = pmin(delta,pnorm);
    		}
    		//           evaluate the function at x + p and calculate its norm.

    		FUNCMULT_EVAL(funcmult,wa2,M,N,wa4);
    		*nfev = *nfev + 1;
    		fnorm1 = enorm(wa4,M);

    		//           compute the scaled actual reduction.

            actred = -one;
            if (p1*fnorm1 < fnorm) {
            	actred = one - (fnorm1/fnorm)*(fnorm1/fnorm);
            }

            //           compute the scaled predicted reduction and
            //           the scaled directional derivative.

            for(j = 0; j < N;++j) {
                wa3[j] = zero;
                l = ipvt[j];
                temp = wa1[l];
                for(i = 0;i <= j;++i) {
                	wa3[i] = wa3[i] + fjac[i*N+j]*temp;
                }
            }

            temp1 = enorm(wa3,N);
            temp1 = temp1/fnorm;
            temp2 = (sqrt(par)*pnorm)/fnorm;
            prered = temp1*temp1 + temp2*temp2/p5;
            dirder = -(temp1*temp1 + temp2*temp2);
            //           compute the ratio of the actual to the predicted
            //           reduction.
            ratio = zero;
            if (prered != zero) {
            	ratio = actred/prered;
            }
            //           update the step bound.
			

            if (ratio <= p25) {//240
            	if (actred >= zero) {
            		temp = p5;
            	}
            	if (actred < zero) {
            		temp = p5*dirder/(dirder + p5*actred);
            	}
            	if (p1*fnorm1 >= fnorm || temp < p1) {
            		temp = p1;
            	}
                delta = temp*pmin(delta,pnorm/p1);
                par = par/temp;
            } else if (par == zero || ratio >= p75){//240 - 260
                delta = pnorm/p5;
                par = p5*par;
            }//260

            //           test for successful iteration.

            if (ratio >= p0001) {//290
            	//           successful iteration. update x, fvec, and their norms.
            	for(j = 0; j < N;++j) {
                    x[j] = wa2[j];
                    wa2[j] = diag[j]*x[j];
            	}
            	for(i = 0; i < M;++i) {
            		fvec[i] = wa4[i];
            	}
                xnorm = enorm(wa2,N);
                fnorm = fnorm1;
                iter = iter + 1;
            }//290
            //           tests for convergence.
            if ((fabs(actred) <= ftol) && (prered <= ftol) && (p5*ratio <= one)) {
            	info = 1;
            }
            if (delta <= xtol*xnorm) {
            	info = 2;
            }
            if ((fabs(actred) <= ftol) && (prered <= ftol) && (p5*ratio <= one) && (info == 2)) {
            	info = 3;
            }
            if (info != 0) {
            	break;
            }

            //           tests for termination and stringent tolerances.
            if (*nfev >= maxfev) {
            	info = 5;
            }
            if ((fabs(actred) <= epsmch) && (prered <= epsmch) && (p5*ratio <= one)) {
            	info = 6;
            }
            if (delta <= epsmch*xnorm) {
            	info = 7;
            }
            if (gnorm <= epsmch) {
            	info = 8;
            }
            if (info != 0) {
            	break;
            }

    	}

        if (info != 0) {
        	break;
        }


    }


    free(wa1);
    free(wa2);
    free(wa3);
    free(wa4);

	return info;
}





