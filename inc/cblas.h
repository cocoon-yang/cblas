/*
 * cblas.h
 *
 *  Created on: 2018-12-24
 *      Author: Chunfeng Yang 
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

enum class MATRIX_TYPE {
	GENERAL = 1,
	SQUARE,
	SYMMETRY,
	TRI_UP,
	TRI_LOW,
};

class MData;


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

/**
@brief   matrix vector multiply
*/
extern "C" MATHLIBRARY_API int cblas_dgemv(char TRANS, int M, int N, double ALPHA, double* A, int LDA,
	double* X, int INCX, double BETA, double* Y, int INCY);

/**
@brief  triangular matrix vector multiply
*/
extern "C" MATHLIBRARY_API void cblas_dtrmv(char uplo, char trans, char diag, int n, double a[], int lda,
	double x[], int incx);

/**
@brief  solving triangular matrix problems
*/
extern "C" MATHLIBRARY_API void cblas_dtrsv(char UPLO, char TRANS, char DIAG, const int  N, double* A,
	const int LDA, double* X, const int INCX);

/**
@brief  performs the rank 1 operation A := alpha*x*y' + A
*/
extern "C" MATHLIBRARY_API int cblas_dger(int M, int N, double ALPHA, double* X, int INCX, double* Y,
		int INCY, double* A, int LDA);



//==============================================================
// Level 3

extern "C" MATHLIBRARY_API void cblas_dgemm(char transa, char transb, int m, int n, int k,
	double alpha, double a[], int lda, double b[], int ldb, double beta,
	double c[], int ldc);

extern "C" MATHLIBRARY_API void cblas_dtrmm(char side, char uplo, char transa, char diag, int m, int n,
	double alpha, double a[], int lda, double b[], int ldb);

extern "C" MATHLIBRARY_API void cblas_dtrsm(char side, char uplo, char transa, char diag, int m, int n,
	double alpha, double a[], int lda, double b[], int ldb);

extern "C" MATHLIBRARY_API void cblas_dpotrf2(char uplo, int n, double *A, int lda, int info);

extern "C" MATHLIBRARY_API void cblas_dgetrf(int m, int n, double* pA, int lda, int* ipiv, int info);

extern "C" MATHLIBRARY_API void cblas_dsyrk(char uplo, char trans, int n, int k,
	double alpha, double* A, int lda, double beta, double* c, int ldc);

//==============================================================
// Lapack

extern "C" MATHLIBRARY_API void cblas_dlas2(double f, double g, double h, double& ssmin, double& ssMAX);

//==============================================================
// Assistant functions
extern "C" MATHLIBRARY_API bool cblas_lsame(char CA, char CB);

extern "C" MATHLIBRARY_API void cblas_dlarf(char side, int m, int n, int l, double* v, int incv, double tau, double* pC, int ldc, double* work);

extern "C" MATHLIBRARY_API void cblas_dlarfg(int N, double& ALPHA, double* X, int INCX, double& TAU);

extern "C" MATHLIBRARY_API void cblas_xerbla(const char* SRNAME, int INFO);

void xerbla(const char* SRNAME, int INFO);

extern "C" MATHLIBRARY_API int ilaenv(int ispec, const char* name, char OPTS, int N1, int N2, int N3, int N4);

extern "C" MATHLIBRARY_API void showMatrix_d(double* data, int rows, int columns);

extern "C" MATHLIBRARY_API double sign(double A, double B);

extern "C" MATHLIBRARY_API double dlamch(const char* CMACH);
double dlamc3(double A, double B);

extern "C" MATHLIBRARY_API void dlasv2(double f, double g, double h, double& ssmin, double& ssmax,
	double& snl, double& csl, double& snr, double& csr);

extern "C" MATHLIBRARY_API double dlapy2(double x, double y); 

extern "C" MATHLIBRARY_API double idamax(int n, double* dx, int incx); 

extern "C" MATHLIBRARY_API int iladlc(int m, int n, double* pA, int lda);

extern "C" MATHLIBRARY_API int iladlr(int m, int n, double* pA, int lda);

#endif /* D_INC_CBLAS_H_ */
