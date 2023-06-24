/*
 * cblas.h
 *
 *  Created on: 2018-12-24
 *      Author: Yang
 */

#ifndef D_INC_CBLAS_H_
#define D_INC_CBLAS_H_
#include <stdbool.h>
#include <stdlib.h>
#include "math.h"
#include <stdio.h>

// #pragma once

#ifdef MATHLIBRARY_EXPORTS
#define MATHLIBRARY_API __declspec(dllexport)
#else
#define MATHLIBRARY_API __declspec(dllimport)
#endif

//==============================================================
// Level 1

extern "C" MATHLIBRARY_API double cblas_dasum(const int N, const double *X, const int incX);

extern "C" MATHLIBRARY_API void cblas_drot(const int N, double *DX, const int INCX, double *DY,
		const int INCY, const double C, const double S);

extern "C" MATHLIBRARY_API void cblas_drotg( double* DA,  double* DB,   double C,  double S);

extern "C" MATHLIBRARY_API void cblas_drotm(const int N, double *DX, const int INCX, double *DY,
		const int INCY, const double* DPARAM);

extern "C" MATHLIBRARY_API void cblas_daxpy(const int N, const double alpha, const double *X,
		const int incX, double *Y, const int incY);

extern "C" MATHLIBRARY_API void cblas_dcopy(const int N, const double *X, const int incX, double *Y,
		const int incY);

extern "C" MATHLIBRARY_API double cblas_ddot(const int N, const double *X, const int incX,
		const double *Y, const int incY);

extern "C" MATHLIBRARY_API double cblas_dnrm2(const int N, const double *X, const int incX);

extern "C" MATHLIBRARY_API void cblas_dscal(const int N, const double DA, double *DX, const int INCX);

extern "C" MATHLIBRARY_API void cblas_dswap(const int N, double *X, const int incX, double *Y,
	const int incY);

extern "C" MATHLIBRARY_API void cblas_dlaswp(int n, double* a, int lda, int k1, int k2, int* ipiv, int incx);

//==============================================================
// Level 2

extern "C" MATHLIBRARY_API int cblas_dger(int M, int N, double ALPHA, double* X, int INCX, double* Y,
		int INCY, double* A, int LDA);

extern "C" MATHLIBRARY_API int cblas_dgemv(char TRANS, int M, int N, double ALPHA, double* A, int LDA,
		double* X, int INCX, double BETA, double* Y, int INCY);

extern "C" MATHLIBRARY_API void cblas_dtrmv(char uplo, char trans, char diag, int n, double a[], int lda,
	double x[], int incx);

//extern "C" MATHLIBRARY_API void cblas_dtrsv(char UPLO, char TRANS, char DIAG, const int  N, double *A, 
//	const int LDA, double *X, const int INCX)
//

//==============================================================
// Level 3

extern "C" MATHLIBRARY_API void cblas_dgemm(char transa, char transb, int m, int n, int k,
	double alpha, double a[], int lda, double b[], int ldb, double beta,
	double c[], int ldc);

extern "C" MATHLIBRARY_API void cblas_dtrmm(char side, char uplo, char transa, char diag, int m, int n,
	double alpha, double a[], int lda, double b[], int ldb);

extern "C" MATHLIBRARY_API void cblas_dtrsm(char side, char uplo, char transa, char diag, int m, int n,
	double alpha, double a[], int lda, double b[], int ldb);


//==============================================================
// Lapack

extern "C" MATHLIBRARY_API void cblas_dlas2(double f, double g, double h, double& ssmin, double& ssMAX);

//==============================================================
// Assistant functions
extern bool cblas_lsame(char CA, char CB);


extern "C" MATHLIBRARY_API void cblas_dlarfg(int N, double& ALPHA, double* X, int INCX, double& TAU);

extern "C" MATHLIBRARY_API void cblas_xerbla(char* SRNAME, int INFO);
void xerbla(char* SRNAME, int INFO);
//int MAX(const int first, const int second );
//double MAX(const double first, const double second);
// int MIN(int first, int second);
#define min(x,y) (((x) < (y)) ? (x) : (y))

#define MAX(x,y) (((x) > (y)) ? (x) : (y))

extern "C" MATHLIBRARY_API void showMatrix_d(double* data, int rows, int columns);

extern "C" MATHLIBRARY_API double sign(double A, double B);


extern "C" MATHLIBRARY_API double dlamch(char* CMACH);
double dlamc3(double A, double B);
extern "C" MATHLIBRARY_API void dlasv2(double f, double g, double h, double& ssmin, double& ssmax,
	double& snl, double& csl, double& snr, double& csr);
extern "C" MATHLIBRARY_API double dlapy2(double x, double y);
#endif /* D_INC_CBLAS_H_ */
