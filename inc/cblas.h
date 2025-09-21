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
#include "mdata.h"
// #pragma once

#ifdef MATHLIBRARY_EXPORTS
#define MATHLIBRARY_API __declspec(dllexport)
#else
#define MATHLIBRARY_API __declspec(dllimport)
#endif

#define C_VEC_OFFSET 1

//enum class MATRIX_TYPE {
//	GENERAL = 1,
//	SQUARE,
//	SYMMETRY,
//	TRI_UP,
//	TRI_LOW,
//};
//
//class MData;


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

/***
  \brief  DSCAL scales a vector by a constant.
  @param[in]  N: int, number of elements in input vector.
  @param[in]  DA: double, specifies the scalar alpha.
  @param[in]  DX: double*, array, dimension [ 1 + ( N - 1 )*abs( INCX ) ] 
  @param[in]  INCX: int, storage spacing between elements of DX.
*/  
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

/**
 * @brief Cholesky factorization 
 * @param[in] uplo char 
 *  = 'U':  Upper triangle of A is stored;
 *  = 'L':  Lower triangle of A is stored.
 * @param[in] n: int,  The order of the matrix A.  N >= 0. 
 * @param A: double*, dimension (LDA,N)  
 *    On entry, the symmetric matrix A.  If UPLO = 'U', the leading
 *    N-by-N upper triangular part of A contains the upper
 *    triangular part of the matrix A, and the strictly lower
 *    triangular part of A is not referenced.  If UPLO = 'L', the
 *    leading N-by-N lower triangular part of A contains the lower
 *    triangular part of the matrix A, and the strictly upper
 *    triangular part of A is not referenced.
 *
 *    On exit, if INFO = 0, the factor U or L from the Cholesky
 *     factorization A = U**T*U or A = L*L**T.
 * @param[in] lda: int, The leading dimension of the array A.  LDA >= max(1,N). 
 * @param[out] info: int, 
 *   = 0:  successful exit
 *   < 0:  if INFO = -i, the i-th argument had an illegal value
 *   > 0:  if INFO = i, the leading principal minor of order i
 *     is not positive, and the factorization could not be
 *     completed.
 * @return 
*/
extern "C" MATHLIBRARY_API void cblas_dpotrf(char uplo, int n, double* A, int lda, int info);

extern "C" MATHLIBRARY_API void cblas_dpotrf2(char uplo, int n, double *A, int lda, int info);

extern "C" MATHLIBRARY_API void cblas_dsyrk(char uplo, char trans, int n, int k,
	double alpha, double* A, int lda, double beta, double* c, int ldc);

//extern "C" MATHLIBRARY_API void dgetrf2(int m, int n, double* pA, int LDA, int* ipiv, int INFO);


//==============================================================
// Lapack



extern "C" MATHLIBRARY_API void clapack_dgesv(int n, int nrhs, double* a, int lda, int* ipiv, double* b,
	int ldb, int info);


/**
 * @brief computes an LU factorization of a general M-by-N matrix A
 *   using partial pivoting with row interchanges.
 * @param[in] m: int, The number of rows of the matrix A.  M >= 0.
 * @param n: int,
 * @param pA: double*,
 * @param lda: int,
 * @param ipiv: int*,
 * @param info: int,
 * @return
*/
extern "C" MATHLIBRARY_API void clapack_dgetrf(int m, int n, double* pA, int lda, int* ipiv, int& info);


extern "C" MATHLIBRARY_API void clapack_dgetrs(char trans, int n, int nrhs, double* a, int lda, int* ipiv,
	double* b, int ldb, int info);


extern "C" MATHLIBRARY_API void cblas_dlas2(double f, double g, double h, double& ssmin, double& ssMAX);

extern "C" MATHLIBRARY_API void clapack_dlarf1f(char side, int m, int n, double* v, int incv, double tau, double* c, int ldc, double* work);

extern "C" MATHLIBRARY_API void clapack_dlarfg(int n, double& alpha, double* x, int incx, double& tau);

extern "C" MATHLIBRARY_API void clapack_dlarf(char side, int m, int n, double*	v, int incv, double tau, double *c, int ldc, double *work);

extern "C" MATHLIBRARY_API void clapack_dgeqr2(int m, int n, double* a, int lda, double* tau, double* work, int& info);

//==============================================================
// Assistant functions
extern "C" MATHLIBRARY_API bool cblas_lsame(char CA, char CB);

extern "C" MATHLIBRARY_API void cblas_dlarf(char side, int m, int n, int l, double* v, int incv, double tau, double* pC, int ldc, double* work);

extern "C" MATHLIBRARY_API void cblas_dlarfg(int N, double& ALPHA, double* X, int INCX, double& TAU);

extern "C" MATHLIBRARY_API void cblas_xerbla(const char* SRNAME, int INFO);

void xerbla(const char* SRNAME, int INFO);

extern "C" MATHLIBRARY_API int ilaenv(int ispec, const char* name, char OPTS, int N1, int N2, int N3, int N4);


extern "C" MATHLIBRARY_API double sign(double A, double B);

extern "C" MATHLIBRARY_API double dlamch(const char* CMACH);
double dlamc3(double A, double B);

extern "C" MATHLIBRARY_API void dlasv2(double f, double g, double h, double& ssmin, double& ssmax,
	double& snl, double& csl, double& snr, double& csr);

extern "C" MATHLIBRARY_API double dlapy2(double x, double y); 

extern "C" MATHLIBRARY_API double cblas_idamax(int n, double* dx, int incx);

extern "C" MATHLIBRARY_API int iladlc(int m, int n, double* pA, int lda);

extern "C" MATHLIBRARY_API int iladlr(int m, int n, double* pA, int lda);

extern "C" MATHLIBRARY_API int cblas_iladlr(int m, int n, double* pA, int lda);

extern "C" MATHLIBRARY_API void showVector_d(double* data, int m);

extern "C" MATHLIBRARY_API void showMatrix_d(double* data, int rows, int columns);

#endif /* D_INC_CBLAS_H_ */
