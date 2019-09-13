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


//==============================================================
// Level 2

extern "C" MATHLIBRARY_API int cblas_dger(int M, int N, double ALPHA, double* X, int INCX, double* Y,
		int INCY, double* A, int LDA);

extern "C" MATHLIBRARY_API int cblas_dgemv(char TRANS, int M, int N, double ALPHA, double* A, int LDA,
		double* X, int INCX, double BETA, double* Y, int INCY);


//==============================================================
// Level 3



//==============================================================
// Assistant functions
bool cblas_lsame(char CA, char CB);
void cblas_xerbla(char* SRNAME, int INFO);
void xerbla(char* SRNAME, int INFO);
int MAX( int first, int second );

#endif /* D_INC_CBLAS_H_ */
