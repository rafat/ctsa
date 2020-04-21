/*
 * matrix.h
 *
 *  Created on: Jul 1, 2013
 *      Author: USER
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define CUTOFF 192
#define TOL 1e-12
#define BLOCKSIZE 64
#define TBLOCK 64
#define SVDMAXITER 50
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define RSVD_POWER_ITERATIONS 5


#ifdef __cplusplus
extern "C" {
#endif

double macheps();

double pmax(double a, double b);

double pmin(double a, double b);

int imax(int a, int b);

int imin(int a, int b);

double signx(double x);

double l2norm(double *vec, int N);

int compare (const void* ind1, const void* ind2);

void sort1d(double* v,int N, int* pos);

//Array Parallel Implementation may have a lot of overhead

double array_max_abs(double *array,int N);

double array_max(double *array,int N);

double array_min(double *array,int N);

//void mmult(double* A, double *B, double *C,int ra,int ca, int rb, int cb);

void dtranspose(double *sig, int rows, int cols,double *col);

void stranspose(double *sig, int rows, int cols,double *col);

void rtranspose(double *m, int rows, int cols,double *n, int r, int c);

void ctranspose(double *sig, int rows, int cols,double *col);

void mtranspose(double *sig, int rows, int cols,double *col);

void itranspose(double *A, int M, int N);

//int minverse(double *xxt, int p);

void mdisplay(double *A, int row, int col);

void madd(double* A, double* B, double* C,int rows,int cols);

void msub(double* A, double* B, double* C,int rows,int cols);

void scale(double *A, int rows, int cols, double alpha);

void nmult(double* A, double* B, double* C,int m,int n, int p);

void tmult(double* A, double* B, double* C,int m,int n, int p);

void recmult(double* A, double* B, double* C,int m,int n, int p,int sA,int sB, int sC);

void rmult(double* A, double* B, double* C,int m,int n, int p);

int findrec(int *a, int *b, int *c);

double house_2(double*x,int N,double *v);

void add_zero_pad(double *X, int rows, int cols, int zrow, int zcol,double *Y);

void remove_zero_pad(double *X, int rows, int cols, int zrow, int zcol,double *Y);

void madd_stride(double* A, double* B, double* C,int rows,int cols,int sA,int sB,int sC);

void msub_stride(double* A, double* B, double* C,int rows,int cols,int sA,int sB,int sC);

void rmadd_stride(double* A, double* B, double* C,int rows,int cols,int p,int sA,int sB,int sC);

void rmsub_stride(double* A, double* B, double* C,int rows,int cols,int p,int sA,int sB,int sC);

void srecmult(double* A, double* B, double* C,int m,int n, int p,int sA,int sB,int sC);

void smult(double* A, double* B, double* C,int m,int n, int p);

void mmult(double* A, double* B, double* C,int m,int n, int p);

void ludecomp(double *A,int N,int *ipiv);

int rludecomp(double *A, int M, int N, int *ipiv);

void getPLU(double *A, int M , int N, int *ipiv,double *P, double *L, double *U);

void getPU(double *A, int M, int N, int *ipiv, double *P,double *U);

double* marsaglia_generate(double *values, int N, double average, double  deviation);

void random_matrix(double *A, int M, int N);

void linsolve(double *A,int N,double *b,int *ipiv,double *x);

void minverse(double *A,int M,int *ipiv,double *inv);

void eye(double *mat,int N);

void eye_scale(double *mat, int N, double lambda);

double house(double*x,int N,double *v);

void housemat(double *v, int N,double beta,double *mat);

void qrdecomp(double *A, int M, int N,double *bvec);

void getQR(double *A,int M,int N,double *bvec,double *Q, double *R);

void hessenberg(double *A,int N);

void francisQR(double *A,int N);

void eig22(double *A, int stride,double *eigre,double *eigim);

int francis_iter(double *A, int N, double *H);

void eig(double *A,int N,double *eigre,double *eigim);

void eigensystem(double *mat, int N, double *eval, double *evec);;

int cholu(double *A, int N);

int bcholu(double *A, int N);

int chol(double *A, int N);

void chold(double *A, int N);

void svd_sort(double *U,int M,int N,double *V,double *q);

int svd(double *A,int M,int N,double *U,double *V,double *q);

int svd_transpose(double *A, int M, int N, double *U, double *V, double *q);

int rank(double *A, int M,int N);

void rsvd(double *A, int M, int N,int K, int oversample, int n_iter,double *U, double *V, double *S);

#ifdef __cplusplus
}
#endif

#endif /* MATRIX_H_ */
