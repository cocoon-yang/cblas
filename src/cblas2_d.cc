#include "cblas.h"
#include <utility>
#include "math.h" 
#include "mdata.h"

extern void showMatrix_d(double* data, int m, int k);
extern void cblas_xerbla(const char* SRNAME, int INFO);
extern bool cblas_lsame(char CA, char CB);
extern void cblas_dlaswp(int n, double* pA, int lda, int k1, int k2, int* ipiv, int incx); 
extern void cblas_dtrsm(char side, char uplo, char transa, char diag, int m, int n,
	double alpha, double* pA, int lda, double* pB, int ldb);
extern void cblas_dgetrf(int m, int n, double* pA, int lda, int* ipiv, int& info); 
//==============================================================
// Level 2

/*
Purpose:

SGER computes A := alpha*x*y' + A.

Discussion:

DGER performs the rank 1 operation

A := alpha*x*y' + A,

where alpha is a scalar, x is an m element vector, y is an n element
vector and A is an m by n matrix.

Licensing:

This code is distributed under the GNU LGPL license.

Modified:

2017

Author:

Original FORTRAN77 version by Jack Dongarra,  Jeremy Du Croz,
Sven Hammarling,  Richard Hanson.
This C version by Chunfeng Yang, John Burkardt.

Parameters:

Input, int M, the number of rows of the matrix A.
0 <= M.

Input, int N, the number of columns of the matrix A.
0 <= N.

Input, float ALPHA, the scalar multiplier.

Input, float X[1+(M-1)*abs(INCX)], the first vector.

Input, int INCX, the increment for elements of X.
INCX must not be zero.

Input, float Y[1+(N-1)*abs(INCY)], the second vector.

Input, int INCY, the increment for elements of Y.
INCY must not be zero.

Input/output, float A[LDA*N].  On entry, the leading M by N
part of the array contains the matrix of coefficients. On exit, A is
overwritten by the updated matrix.

Input, int LDA, the first dimension of A as declared
in the calling program. max ( 1, M ) <= LDA.
*/
// row prinary order version

/**
  @brief      GER computes A := alpha*x*y' + A. 
  @warning    This is ROW prinary order version 
  @param[in]  M: int,   the number of rows of the matrix A. 0 <= M. 
  @param[in]  N: int,   the number of columns of the matrix A. 0 <= M.
 ***/
int cblas_dger(int M, int N, double ALPHA, double* X, int INCX, double* Y,
	int INCY, double* A, int LDA)
{
	// local variables
	double TEMP;
	int I, INFO, IX, J, JY, KX;

	double ZERO = 0.0;

	/*
	*     Test the input parameters.
	*/
	INFO = 0;
	if (M < 0)
	{
		INFO = 1;
	}
	else if (N < 0)
	{
		INFO = 2;
	}
	else if (INCX == 0)
	{
		INFO = 5;
	}
	else if (INCY == 0)
	{
		INFO = 7;
	}
	else if (LDA < std::max(1, N))
	{
		INFO = 9;
	}
	if (INFO != 0)
	{
		xerbla("DGER ", INFO);
		return 0;
	}
	/*
	*     Quick return if possible.
	*/
	if ((M == 0) || (N == 0) || (ALPHA == ZERO))
	{
		return 0;
	}

	/*
	*     Start the operations. In this version the elements of A are
	*     accessed sequentially with one pass through A.
	*/
	if (INCY > 0)
	{
		JY = 0;
	}
	else
	{
		JY = 0 - (N - 1) * INCY;
	}

	// 
	// Row Prinary Order -- BEGIN --
	if (INCX == 1)
	{
		for (I = 0; I < M; I++)
		{
			if (X[I] != ZERO)
			{
				JY = 0;
				TEMP = ALPHA * X[I];
				for (J = 0; J < N; J = J + INCY)
				{
					A[(I)*LDA + (J)] = A[(I)*LDA + (J)] + Y[JY] * TEMP;
					JY += INCY;
				}
			}
			//JY = JY + INCY;
		}
	}
	else
	{
		IX = 0;
		for (I = 0; I < M; I++)
		{
			if (X[IX] != ZERO)
			{
				JY = 0;
				TEMP = ALPHA * X[IX];
				for (J = 0; J < N; J++)
				{
					A[(I)*LDA + (J)] = A[(I)*LDA + (J)] + Y[JY] * TEMP;
					JY += INCY;
				}
			}
			IX += INCX;
		}
	}
	// Row Prinary Order -- END -- 
	//
	return 0;
}



/**
@brief      GEMV computes A := alpha*A*x + A.
@warning    This is row prinary order version
@param[in]  M: int,   the number of rows of the matrix A. 0 <= M.
@param[in]  N: int,   the number of columns of the matrix A. 0 <= M.
***/
int cblas_dgemv(char TRANS, int M, int N, double ALPHA, double* pA, int LDA,
	double* X, int INCX, double BETA, double* Y, int INCY)
{
	double ONE = 1.0;
	double ZERO = 0.0;

	double TEMP;
	int I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY;

	INFO = 0;

	if (!cblas_lsame(TRANS, 'N') && !cblas_lsame(TRANS, 'T') && !cblas_lsame(
		TRANS, 'C'))
	{
		INFO = 1;
	}
	else if (M < 0)
	{
		INFO = 2;
	}
	else if (N < 0)
	{
		INFO = 3;
	}
	else if (LDA < std::max(1, N))
	{
		INFO = 6;
	}
	else if (INCX == 0)
	{
		INFO = 8;
	}
	else if (INCY == 0)
	{
		INFO = 11;
	}
	if (INFO != 0)
	{
		xerbla("DGEMV ", INFO);
		return 0;
	}

	/*
	*     Quick return if possible.
	*/
	//	      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
	//	     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
	if ((M == 0) || (N == 0) || ((ALPHA == ZERO) && (BETA == ONE)))
	{
		return 0;
	}

	/*
	* Set  LENX  and  LENY, the lengths of the vectors x and y, and set
	*  up the start points in  X  and  Y.
	*/

	if (cblas_lsame(TRANS, 'N'))
	{
		LENX = N;
		LENY = M;
	}
	else
	{
		LENX = M;
		LENY = N;
	}

	if (INCX > 0)
	{
		//KX = 1;  // Fortran version
		KX = 0; // C/C++ version
	}
	else
	{
		//KX = 1 - (LENX-1)*INCX;// Fortran version
		KX = -(LENX - 1) * INCX; // C/C++ version
	}

	if (INCY > 0)
	{
		//KY = 1; // Fortran version
		KY = 0; // C/C++ version
	}
	else
	{
		//KY = 1 - (LENY-1)*INCY; // Fortran version
		KY = -(LENY - 1) * INCY; // C/C++ version
	}

	/*
	*     Start the operations. In this version the elements of A are
	*     accessed sequentially with one pass through A.
	*
	*     First form  y := beta*y.
	*/
	if (BETA != ONE)
	{
		if (INCY == 1)
		{
			// if parameter Beta is 0
			// set vector Y = [0]
			if (BETA == ZERO)
			{
				// for (I = 1; I < LENY; I++)  // Fortran version
				for (I = 0; I < LENY; I++) // C/C++ version
				{
					Y[I] = ZERO;
				}
			}

			else
			{
				// for (I = 1; I < LENY; I++)  // Fortran version
				for (I = 0; I < LENY; I++) // C/C++ version
				{
					Y[I] = BETA * Y[I];
				}
			}
		}
		else // (INCY != 1)
		{
			IY = KY;
			if (BETA == ZERO)
			{
				// for (I = 1; I < LENY; I++)  // Fortran version
				for (I = 0; I < LENY; I++) // C/C++ version
				{
					Y[IY] = ZERO;
					IY = IY + INCY;
				}
			}
			else
			{
				// for (I = 1; I < LENY; I++)  // Fortran version
				for (I = 0; I < LENY; I++) // C/C++ version
				{
					Y[IY] = BETA * Y[IY];
					IY = IY + INCY;
				}
			}
		}
	}

	//    IF (ALPHA.EQ.ZERO) RETURN
	if (ALPHA == ZERO)
	{
		return 0;
	}

	MData A(M, LDA, pA, LDA);
	//A.setData(pA);

	if (cblas_lsame(TRANS, 'N'))
	{
		// 
		//  Form  y := alpha*A*x + y.
		// 
		JX = KX;
		if (INCY == 1)
		{
			for (J = 0; J < N; J++)
			{
				if (X[JX] != ZERO)
				{
					TEMP = ALPHA * X[JX];
					for (I = 0; I < M; I++)
					{
						Y[I] = Y[I] + TEMP * A[I][J];
					}
				}
				JX = JX + INCX;
			}
		}
		else
		{
			for (J = 0; J < N; J++)
			{
				if (X[JX] != ZERO)
				{
					TEMP = ALPHA * X[JX];
					IY = KY;
					for (I = 0; I < M; I++)
						// for (I = 1; I < M; I++) // Fortran version
					{
						Y[IY] = Y[IY] + TEMP * A[I][J];
						IY = IY + INCY;
					}
				}
				JX = JX + INCX;
			}
		}
	}
	else
	{
		/*
		*        Form  y := alpha*A'*x + y.
		*/
		JY = KY;
		if (INCX == 1)
		{
			for (J = 0; J < N; J++)
			{
				TEMP = ZERO;
				for (I = 0; I < M; I++)
				{
					TEMP = TEMP + A[I][J] * X[I];
					// std::cout << "TEMP = TEMP + A[" << I << "][" << J << "] * X[" << I << "]" << std::endl;
				}
				Y[JY] = Y[JY] + ALPHA * TEMP;
				JY = JY + INCY;
			}
		}
		else
		{
			for (J = 0; J < N; J++)
			{
				TEMP = ZERO;
				IX = KX;
				for (I = 0; I < M; I++)
				{
					TEMP = TEMP + A[I][J] * X[IX];
					IX = IX + INCX;
				}
				Y[JY] = Y[JY] + ALPHA * TEMP;
				JY = JY + INCY;
			}
		}
	}

}





//****************************************************************************80
//
//  Purpose:
//
//    DTRMV computes x: = A*x or x = A'*x for a triangular matrix A.
//
//  Discussion:
//
//    DTRMV performs one of the matrix-vector operations
//
//      x := A*x,   or   x := A'*x,
//
//    where x is an n element vector and  A is an n by n unit, or non-unit,
//    upper or lower triangular matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2014
//
//  Author:
//
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, char UPLO, specifies whether the matrix is an upper or
//    lower triangular matrix as follows:
//    'u' or 'U': A is an upper triangular matrix.
//    'l' or 'L': A is a lower triangular matrix.
//
//    Input, char TRANS, specifies the operation to be performed as
//    follows:
//    'n' or 'N': x := A*x.
//    't' or 'T': x := A'*x.
//    'c' or 'C': x := A'*x.
//
//    Input, char DIAG, specifies whether or not A is unit
//    triangular as follows:
//    'u' or 'U': A is assumed to be unit triangular.
//    'n' or 'N': A is not assumed to be unit triangular.
//
//    Input, int N, the order of the matrix A.
//    0 <= N.
//
//    Input, double A[LDA*N].
//    Before entry with  UPLO = 'u' or 'U', the leading n by n
//    upper triangular part of the array A must contain the upper
//    triangular matrix and the strictly lower triangular part of
//    A is not referenced.
//    Before entry with UPLO = 'l' or 'L', the leading n by n
//    lower triangular part of the array A must contain the lower
//    triangular matrix and the strictly upper triangular part of
//    A is not referenced.
//    Note that when  DIAG = 'u' or 'U', the diagonal elements of
//    A are not referenced either, but are assumed to be unity.
//
//    Input, int LDA, the first dimension of A as declared
//    in the calling program. max ( 1, N ) <= LDA.
//
//    Input/output, double X[1+(N-1)*abs( INCX)].
//    Before entry, the incremented array X must contain the n
//    element vector x. On exit, X is overwritten with the
//    tranformed vector x.
//
//    Input, int INCX, the increment for the elements of
//    X.  INCX must not be zero.
//
void cblas_dtrmv(char uplo, char trans, char diag, int n, double a[], int lda,
	double x[], int incx)
{
	int i;
	int info;
	int ix;
	int j;
	int jx;
	int kx;
	bool nounit;
	double temp;
	//
	//  Test the input parameters.
	//
	info = 0;
	if (!cblas_lsame(uplo, 'U') && !cblas_lsame(uplo, 'L'))
	{
		info = 1;
	}
	else if (!cblas_lsame(trans, 'N') && !cblas_lsame(trans, 'T') &&
		!cblas_lsame(trans, 'C'))
	{
		info = 2;
	}
	else if (!cblas_lsame(diag, 'U') && !cblas_lsame(diag, 'N'))
	{
		info = 3;
	}
	else if (n < 0)
	{
		info = 4;
	}
	else if (lda < std::max(1, n))
	{
		info = 6;
	}
	else if (incx == 0)
	{
		info = 8;
	}

	if (info != 0)
	{
		//xerbla("DTRMV", info);
		return;
	}
	//
	//  Quick return if possible.
	//
	if (n == 0)
	{
		return;
	}

	nounit = cblas_lsame(diag, 'N');
	//
	//  Set up the start point in X if the increment is not unity. This
	//  will be  ( N - 1 ) * INCX  too small for descending loops.
	//
	if (incx <= 0)
	{
		kx = 0 - (n - 1) * incx;
	}
	else if (incx != 1)
	{
		kx = 0;
	}

	MData A(n, lda, a, lda);
	//A.setData(a);

	//
	//  Start the operations. In this version the elements of A are
	//  accessed sequentially with one pass through A.
	//
	if (cblas_lsame(trans, 'N'))
	{
		//
		//  Form x := A*x.
		//
		if (cblas_lsame(uplo, 'U'))
		{
			if (incx == 1)
			{
				for (j = 0; j < n; j++)
				{
					if (x[j] != 0.0)
					{
						temp = x[j];
						for (i = 0; i < j; i++)
						{
							x[i] = x[i] + temp * A[i][j];
						}
						if (nounit)
						{
							x[j] = x[j] * A[j][j];
						}
					}
				}
			}
			else
			{
				jx = kx;
				for (j = 0; j < n; j++)
				{
					if (x[jx] != 0.0)
					{
						temp = x[jx];
						ix = kx;
						for (i = 0; i < j; i++)
						{
							x[ix] = x[ix] + temp * A[i][j];
							ix = ix + incx;
						}
						if (nounit)
						{
							x[jx] = x[jx] * A[j][j];
						}
					}
					jx = jx + incx;
				}
			}
		}
		else
		{
			if (incx == 1)
			{
				for (j = n - 1; 0 <= j; j--)
				{
					if (x[j] != 0.0)
					{
						temp = x[j];
						for (i = n - 1; j < i; i--)
						{
							x[i] = x[i] + temp * A[i][j];
						}
						if (nounit)
						{
							x[j] = x[j] * A[j][j];
						}
					}
				}
			}
			else
			{
				kx = kx + (n - 1) * incx;
				jx = kx;
				for (j = n - 1; 0 <= j; j--)
				{
					if (x[jx] != 0.0)
					{
						temp = x[jx];
						ix = kx;
						for (i = n - 1; j < i; i--)
						{
							x[ix] = x[ix] + temp * A[i][j];
							ix = ix - incx;
						}
						if (nounit)
						{
							x[jx] = x[jx] * A[j][j];
						}
					}
					jx = jx - incx;
				}
			}
		}
	}
	else
	{
		//
		//  Form x := A'*x.
		//
		if (cblas_lsame(uplo, 'U'))
		{
			if (incx == 1)
			{
				for (j = n - 1; 0 <= j; j--)
				{
					temp = x[j];
					if (nounit)
					{
						temp = temp * A[j][j];
					}
					for (i = j - 1; 0 <= i; i--)
					{
						temp = temp + A[i][j] * x[i];
					}
					x[j] = temp;
				}
			}
			else
			{
				jx = kx + (n - 1) * incx;
				for (j = n - 1; 0 <= j; j--)
				{
					temp = x[jx];
					ix = jx;
					if (nounit)
					{
						temp = temp * A[j][j];
					}
					for (i = j - 1; 0 <= i; i--)
					{
						ix = ix - incx;
						temp = temp + A[i][j] * x[ix];
					}
					x[jx] = temp;
					jx = jx - incx;
				}
			}
		}
		else
		{
			if (incx == 1)
			{
				for (j = 0; j < n; j++)
				{
					temp = x[j];
					if (nounit)
					{
						temp = temp * A[j][j];
					}
					for (i = j + 1; i < n; i++)
					{
						temp = temp + A[i][j] * x[i];
						//std::cout << "temp = temp + A[" << i << "][" << j << "] * x[" << i << "]; " << std::endl;
					}
					x[j] = temp;
				}
			}
			else
			{
				jx = kx;
				for (j = 0; j < n; j++)
				{
					temp = x[jx];
					ix = jx;
					if (nounit)
					{
						temp = temp * A[j][j];
					}
					for (i = j + 1; i < n; i++)
					{
						ix = ix + incx;
						temp = temp + A[i][j] * x[ix];
					}
					x[jx] = temp;
					jx = jx + incx;
				}
			}
		}
	}
	return;
}


/***
@brief
 DTRSV  solves one of the systems of equations

	A*x = b,   or   A**T*x = b,

 where b and x are n element vectors and A is an n by n unit, or
 non-unit, upper or lower triangular matrix.

 No test for singularity or near-singularity is included in this
 routine. Such tests must be performed before calling this routine.

 @param[in] uplo: char,
			On entry, UPLO specifies whether the matrix is an upper or
		   lower triangular matrix as follows:

			  UPLO = 'U' or 'u'   A is an upper triangular matrix.

			  UPLO = 'L' or 'l'   A is a lower triangular matrix.

 @param[in] trans: char,
			On entry, TRANS specifies the equations to be solved as
		   follows:

			  TRANS = 'N' or 'n'   A*x = b.

			  TRANS = 'T' or 't'   A**T*x = b.

			  TRANS = 'C' or 'c'   A**T*x = b.

 @param[in] diag: char,
			On entry, DIAG specifies whether or not A is unit
		   triangular as follows:

			  DIAG = 'U' or 'u'   A is assumed to be unit triangular.

			  DIAG = 'N' or 'n'   A is not assumed to be unit
								  triangular.

 @param[in] n: int,
			On entry, N specifies the order of the matrix A.
		   N must be at least zero.

 @param[in] pA: double*
			Before entry with  UPLO = 'U' or 'u', the leading n by n
		   upper triangular part of the array A must contain the upper
		   triangular matrix and the strictly lower triangular part of
		   A is not referenced.
		   Before entry with UPLO = 'L' or 'l', the leading n by n
		   lower triangular part of the array A must contain the lower
		   triangular matrix and the strictly upper triangular part of
		   A is not referenced.
		   Note that when  DIAG = 'U' or 'u', the diagonal elements of
		   A are not referenced either, but are assumed to be unity.

 @param[in] lda: int
			On entry, LDA specifies the first dimension of A as declared
		   in the calling (sub) program. LDA must be at least
		   max( 1, n ).

 @param[in] x: double*
			Before entry, the incremented array X must contain the n
		   element right-hand side vector b. On exit, X is overwritten
		   with the solution vector x.

 @param[in] incx: int
			On entry, INCX specifies the increment for the elements of
		   X. INCX must not be zero.
*/
void cblas_dtrsv(char uplo, char trans, char diag, int n, double* pA, int lda, double* x, int incx)
{
	double zero = 0.0;

	double temp;
	int i, info, ix, j, jx, kx;

	bool nounit;
	/*
	* Test the input parameters.
	*/

	info = 0;
	if (!cblas_lsame(uplo, 'U') && !cblas_lsame(uplo, 'L')) {
		info = 1;
	}
	else if (!cblas_lsame(trans, 'N') && !cblas_lsame(trans, 'T') && !cblas_lsame(trans, 'C')) {
		info = 2;
	}
	else if (!cblas_lsame(diag, 'U') && !cblas_lsame(diag, 'N')) {
		info = 3;
	}
	else if (n < 0) {
		info = 4;
	}
	else if (lda < std::max(1, n)) {
		info = 6;
	}
	else if (incx == 0) {
		info = 8;
	}
	if (info != 0) {
		cblas_xerbla("DTRSV ", info);
		return;
	}

	/*
	* Quick return if possible.
	*/
	if (n == 0)
	{
		return;
	}

	nounit = cblas_lsame(diag, 'N');

	/*
	*     Set up the start point in X if the increment is not unity. This
	*     will be  ( N - 1 )*INCX  too small for descending loops.
	*/

	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	}
	else if (incx != 1) {
		kx = 1;
	}

	MData a(n, lda, pA, lda);
	//a.setData(pA);


	/*
	*     Start the operations. In this version the elements of A are
	*     accessed sequentially with one pass through A.
	*/

	if (cblas_lsame(trans, 'N')) {
		/*
		*        Form  x := inv( A )*x.
		*/
		if (cblas_lsame(uplo, 'U'))
		{
			if (incx == 1) {
				for (j = n - 1; j >= 0; j--)
				{
					if (x[j] != zero) {
						if (nounit)
						{
							x[j] = x[j] / a[j][j];
						}
						temp = x[j];
						for (i = j - 1; i >= 0; i--)
						{
							x[i] = x[i] - temp * a[i][j];
						}
					}
				}
			}
			else {
				jx = kx + (n - 1) * incx;
				for (j = n - 1; j >= 0; j--)
				{
					if (x[jx] != zero)
					{
						if (nounit)
						{
							x[jx] = x[jx] / a[j][j];
						}
						temp = x[jx];
						ix = jx;
						for (i = j - 1; i >= 0; i--)
						{
							ix = ix - incx;
							x[ix] = x[ix] - temp * a[i][j];
						}
					}
					jx = jx - incx;
				}
			}
		}
		else { // 'L'

			if (incx == 1)
			{
				for (j = 0; j < n; j++)
				{
					if (x[j] != zero)
					{
						if (nounit)
						{
							x[j] = x[j] / a[j][j];
						}
						temp = x[j];
						for (i = j + 1; i < n; i++)
						{
							x[i] = x[i] - temp * a[i][j];
						}
					}
				}
			}
			else {
				jx = kx;
				for (j = 0; j < n; j++)
				{
					if (x[jx] != zero)
					{
						if (nounit)
						{
							x[jx] = x[jx] / a[j][j];
						}
						temp = x[jx];
						ix = jx;
						for (i = j + 1; i < n; i++)
						{
							ix = ix + incx;
							x[ix] = x[ix] - temp * a[i][j];
						}
					}
					jx = jx + incx;
				}
			}
		}
	}
	else {
		/*
		*        Form  x := inv( A**T )*x.
		*/
		if (cblas_lsame(uplo, 'U'))
		{
			if (incx == 1) {
				for (j = 0; j < n; j++)
				{
					temp = x[j];
					for (i = 0; i < j; i++)
					{
						temp = temp - a[i][j] * x[i];
					}
					if (nounit)
					{
						temp = temp / a[j][j];
					}
					x[j] = temp;
				}
			}
			else {
				jx = kx;
				for (j = 0; j < n; j++)
				{
					temp = x[jx];
					ix = kx;
					for (i = 0; i < j; i++)
					{
						temp = temp - a[i][j] * x[ix];
						ix = ix + incx;
					}
					if (nounit)
					{
						temp = temp / a[j][j];
					}
					x[jx] = temp;
					jx = jx + incx;
				}
			} 
		}
		else {
			if (incx == 1) {
				for (j = n - 1; j >= 0; j--)
				{
					temp = x[j];
					for (i = n - 1; i >= (j + 1); i--)
					{
						temp = temp - a[i][j] * x[i];
					}
					if (nounit)
					{
						temp = temp / a[j][j];
					}
					x[j] = temp;
				}
			}
			else {
				kx = kx + (n - 1) * incx;
				jx = kx;
				for (j = n - 1; j >= 0; j--)
				{
					temp = x[jx];
					ix = kx;
					for (i = n - 1; i >= (j + 1); i--)
					{
						temp = temp - a[i][j] * x[ix];
						ix = ix - incx;
					}
					if (nounit)
					{
						temp = temp / a[j][j];
					}
					x[jx] = temp;
					jx = jx - incx;
				}
			}
		}
	}
}

 
//
// DTPMV  performs one of the matrix-vector operations
//
//    x := A*x,   or   x := A**T*x,
//
// where x is an n element vector and  A is an n by n unit, or non-unit,
// upper or lower triangular matrix, supplied in packed form.
//
//          UPLO is CHARACTER*1
//           On entry, UPLO specifies whether the matrix is an upper or
//           lower triangular matrix as follows:
//
//              UPLO = 'U' or 'u'   A is an upper triangular matrix.
//
//              UPLO = 'L' or 'l'   A is a lower triangular matrix.
//
//          AP is DOUBLE PRECISION array, dimension at least
//           ( ( n*( n + 1 ) )/2 ).
//           Before entry with  UPLO = 'U' or 'u', the array AP must
//           contain the upper triangular matrix packed sequentially,
//           column by column, so that AP( 1 ) contains a( 1, 1 ),
//           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
//           respectively, and so on.
// 
//   In Fortran version, when UPLO = 'U', the upper triangle matrix
//  [ a0, a1, a3, a6, ]
//  [     a2, a4, a7, ]
//  [         a5, a8, ]
//  [             a9  ]  
// 
//  
//           Before entry with UPLO = 'L' or 'l', the array AP must
//           contain the lower triangular matrix packed sequentially,
//           column by column, so that AP( 1 ) contains a( 1, 1 ),
//           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
//           respectively, and so on.
// 
// 
//   In Fortran version, when UPLO = 'L', the lower triangle matrix
//  [ a0,            ]
//  [ a1, a4,        ]
//  [ a2, a5, a7,    ]
//  [ a3, a6, a8, a9 ]  
// 
//           Note that when  DIAG = 'U' or 'u', the diagonal elements of
//           A are not referenced, but are assumed to be unity.
//
void cblas_dtpmv(char  uplo, char trans, char diag, int n, double* ap, double* x, int incx)
{
	//
	//     Test the input parameters.
	//
	int type = 4;

	int info = 0;
	if (!cblas_lsame(uplo, 'U') && !cblas_lsame(uplo, 'L')) {
		info = 1;
	}
	else if (!cblas_lsame(trans, 'N') &&
		!cblas_lsame(trans, 'T') &&
		!cblas_lsame(trans, 'C')) {
		info = 2;
	}
	else if (!cblas_lsame(diag, 'U') &&
		!cblas_lsame(diag, 'N')) {
		info = 3;
	}
	else if (n < 0) {
		info = 4;
	}
	else if (incx == 0) {
		info = 7;
	}
	if (info != 0) {
		cblas_xerbla("DTPMV ", info);
		return;
	}
	//
	//     Quick return if possible.
	//
	if (n == 0) return;

	double zero = 0.0;
	double temp;
	int k, kk, ix, jx, kx;

	if (cblas_lsame(uplo, 'U'))
	{
		type = 4;
	}
	else if (cblas_lsame(uplo, 'L'))
	{
		type = 5;
	}
	MData t(n, n, ap, n, type);
	//t.show(); 
	//double foo = t(0,2);

	//
	bool nounit = cblas_lsame(diag, 'N');
	//
	//     Set up the start point in X if the increment is not unity. This
	//     will be  ( N - 1 )*INCX  too small for descending loops.
	//
	if (incx <= 0) {
		kx = -(n - 1) * incx;
	}
	else if (incx != 1) {
		kx = 0;
	}
	//
	//     Start the operations. In this version the elements of AP are
	//     accessed sequentially with one pass through AP.
	//
	if (cblas_lsame(trans, 'N')) {
		//
		//        Form  x:= A*x.
		// 
		if (cblas_lsame(uplo, 'U')) {
			//
			//  [ a0, a1, a2, a3, ]
			//  [     a4, a5, a6, ]
			//  [         a7, a8, ]
			//  [             a9  ]  
			// 
			kk = 1;
			if (incx == 1) {
				k = 0;
				for (int j = 0; j < n; j++)
				{
					temp = x[j];
					if (nounit)
					{
						x[j] = x[j] * ap[k];
						//
						// printf("  x[%i] = x[%i] * ap[%i] = %5.3f \n", j, j, k, x[j]);  
						// 
					}
					k++;

					for (int i = j + 1; i < n; i++)
					{
						x[j] = x[j] + x[i] * ap[k];

						//
						// printf("  x[%i] = x[%i] + x[%i] * ap[%i] = %5.3f \n", j, j, i, k, x[j]);
						//  
						k++;
					}
				}


				//for (int j = 0; j < n; j++) {
				//	if (x[j] != zero) {
				//		temp = x[j];
				//		k = kk;
				//		for (int i = 0; i < j; i++) {
				//			//x[i] = x[i] + temp * ap[k]; 
				//			x[i] = x[i] + temp * t(j,i);

				//			printf("  x[%i] = x[%i] + temp * t[%i][%i]  \n", i, i, i, j);
				//			printf("   t[%i][%i]  = %5.3f \n", i, j, t(j, i));
				//				 

				//			k = k + 1;
				//		}
				//		if (nounit)
				//		{
				//			// x[j] = x[j] * ap[kk + j - 1]; 
				//			x[j] = x[j] * t(j, j); // [j] [j] ;
				//			printf("  x[%i] = x[%i] * ap[%i]  \n", j, j, kk + j - 1);
				//		}
				//	}
				//	kk = kk + j;
				//}
			}
			else {
				k = 0;
				jx = kx;
				for (int j = 0; j < n; j++)
				{
					temp = x[jx];
					if (nounit)
					{
						x[jx] = x[jx] * ap[k];
						//printf("  x[%i] = x[%i] * ap[%i] = %5.3f \n", j, j, k, x[j]);
					}
					k++;

					ix = kx;
					for (int i = j + 1; i < n; i++)
					{
						x[jx] = x[jx] + x[ix] * ap[k];

						//printf("  x[%i] = x[%i] + x[%i] * ap[%i] = %5.3f \n", j, j, i, k, x[j]);
						ix = ix + incx;
						k++;
					}
					jx = jx + incx;

				}


				//jx = kx;
				//for (int j = 0; j < n; j++) {
				//	if (x[jx] != zero) {
				//		temp = x[jx];
				//		ix = kx;
				//		for (int k = kk - 1; k < kk + j - 2; k++) {
				//			x[ix] = x[ix] + temp * ap[k];
				//			ix = ix + incx;
				//		}
				//		if (nounit) x[jx] = x[jx] * ap[kk + j - 1];
				//	}
				//	jx = jx + incx;
				//	kk = kk + j;
				//}
			}


		}
		else {
			// uplo -- 'L'
			//
			//  [ a0,        ]
			//  [ a1, a2,    ]
			//  [ a3, a4, a5 ]
			//  [ a6, ...    ]  
			// 

			kk = (n * (n + 1)) / 2 - 1;
			if (incx == 1) {
				k = kk;
				for (int j = n - 1; j >= 0; j--) {
					if (x[j] != zero) {
						temp = x[j];

						if (nounit)
						{
							x[j] = x[j] * ap[k];
						}
						k = k - 1;
						for (int i = j - 1; i >= 0; i--) {
							x[j] = x[j] + x[i] * ap[k];
							k = k - 1;
						}
					}
				}
			}
			else {
				k = kk;
				kx = kx + (n - 1) * incx;
				jx = kx;
				for (int j = n - 1; j >= 0; j--) {
					if (x[jx] != zero) {
						temp = x[jx];

						if (nounit)
						{
							x[jx] = x[jx] * ap[k];
						}
						k = k - 1;
						for (int i = j - 1; i >= 0; i--) {
							x[jx] = x[jx] + x[i] * ap[k];
							k = k - 1;
						}
					}
					jx = jx - incx;
				}

				//for (int j = n - 1; j >= 0; j--) {
				//	if (x[jx] != zero) {
				//		temp = x[jx];
				//		ix = kx;
				//		for (int k = kk - 1; k >= kk - (n - (j + 1)); k--) {
				//			x[ix] = x[ix] + temp * ap[k];
				//			ix = ix - incx;
				//		}
				//		if (nounit) 
				//		{ 
				//			x[jx] = x[jx] * ap[kk - n + j];
				//		}
				//	}
				//	jx = jx - incx;
				//	kk = kk - (n - j + 1);
				//}
			}
		}
	}
	else {
		//
		//  Form  x := A**T*x.
		//
		if (cblas_lsame(uplo, 'U')) {
			//
			//  Upper Triangle 
			// 
			//  [ a0, a1, a2, a3, ]
			//  [     a4, a5, a6, ]
			//  [         a7, a8, ]
			//  [             a9  ]  
			//  
			//  After transpose, the matrix is
			// 
			//  [ a0,            ]
			//  [ a1, a4,        ]
			//  [ a2, a5, a7,    ]
			//  [ a3, a6, a8, a9 ]  
			//   
			//  x3 = x0*a3 + x1*a6 + x2*a8 + x3*a9 
			//  x2 = x0*a2 + x1*a5 + x2*a7 
			//  x1 = x0*a1 + x1*a4 
			//  x0 = x0*a0 
			//
			//
			kk = (n * (n + 1)) / 2 - 1;
			if (incx == 1) {
				for (int j = n - 1; j >= 0; j--) {
					temp = x[j];
					if (nounit) {
						temp = temp * ap[kk];
					}
					k = kk;
					for (int i = j - 1; i >= 0; i--)
					{
						k = k - (n - i - 1);
						//
						//printf("  temp = temp + ap[%i] * x[%i]", k, i);
						//

						temp = temp + ap[k] * x[i];
						//
						//printf(" = %5.3f \n", temp);
						//  


						//temp = temp + ap[k] * x[i]; 
						//k = k - 1;
					}
					x[j] = temp;
					kk = kk - (n - j + 1);
				}
			}
			else {
				jx = kx + (n - 1) * incx;
				for (int j = n - 1; j >= 0; j--) {
					temp = x[jx];
					ix = jx;
					if (nounit) {
						temp = temp * ap[kk];
					}
					k = kk;
					for (int i = j - 1; i >= 0; i--)
					{
						k = k - (n - i - 1);
						ix = ix - incx;
						temp = temp + ap[k] * x[ix];

					}
					x[jx] = temp;
					kk = kk - (n - j + 1);
					jx = jx - incx;

					//for (int k = kk - 1; k >= kk - j + 1; k--) {
					//	ix = ix - incx;
					//	temp = temp + ap[k] * x[ix];
					//}
					//x[jx] = temp;
					//jx = jx - incx;
					//kk = kk - j;
				}
			}
		}
		else {
			//
			//  Lower Triangle 
			// 
			//  [ a0,         ]
			//  [ a1, a2,     ]
			//  [ a3, a4, a5, ]
			//  [ a6, ...     ]  
			//  
			//  After transpose, the matrix is
			//  [ a0, a1, a3, a6, ]    [x0] 
			//  [     a2, a4, a7, ]    [x1] 
			//  [         a5, a8, ]    [x2] 
			//  [             a9  ]    [x3]   
			// 
			//  x0 = x0*a0 + x1*a1 + x2*a3 + x3*a6
			//  x1 =         x1*a2 + x2*a4 + x3*a7
			//  x2 =                 x2*a5 + x3*a8
			//  x3 =                         x3*a9
			kk = 0;
			if (incx == 1) {
				for (int j = 0; j < n; j++) {
					temp = x[j];
					if (nounit) {
						temp = temp * ap[kk];
					}
					k = kk;
					for (int i = j + 1; i < n; i++) {

						k = k + i;
						//
						//printf("  temp = temp + ap[%i] * x[%i]", k, i);
						//

						temp = temp + ap[k] * x[i];
						//
						// printf(" = %5.3f \n", temp);
						//  


					}
					x[j] = temp;
					kk = kk + (j + 2);
				}
			}
			else {
				jx = kx;
				for (int j = 0; j < n; j++) {
					temp = x[jx];
					ix = jx;
					if (nounit) {
						temp = temp * ap[kk];
					}
					k = kk;
					for (int i = j + 1; i < n; i++)
					{
						ix = ix + incx;
						k = k + i;
						temp = temp + ap[k] * x[ix];
					}
					x[jx] = temp;
					kk = kk + (j + 2);
					jx = jx + incx;

					//for (int k = kk; k < kk + n - j; k++) {
					//	ix = ix + incx;
					//	temp = temp + ap[k] * x[ix];
					//}
					//x[jx] = temp;
					//jx = jx + incx;
					//kk = kk + (n - j + 1);
				}
			}
		}
	}

	return;
}

//
// DTPSV  solves one of the systems of equations
//
//    A*x = b,   or   A**T*x = b,
//
// where b and x are n element vectors and A is an n by n unit, or
// non-unit, upper or lower triangular matrix, supplied in packed form.
//
// No test for singularity or near-singularity is included in this
// routine. Such tests must be performed before calling this routine.
//
//          AP is DOUBLE PRECISION array, dimension at least
//           ( ( n*( n + 1 ) )/2 ).
//           Before entry with  UPLO = 'U' or 'u', the array AP must
//           contain the upper triangular matrix packed sequentially,
//           column by column, so that AP( 1 ) contains a( 1, 1 ),
//           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
//           respectively, and so on.
//           Before entry with UPLO = 'L' or 'l', the array AP must
//           contain the lower triangular matrix packed sequentially,
//           column by column, so that AP( 1 ) contains a( 1, 1 ),
//           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
//           respectively, and so on.
//           Note that when  DIAG = 'U' or 'u', the diagonal elements of
//           A are not referenced, but are assumed to be unity.
//
void cblas_dtpsv(char uplo, char trans, char diag, int n, double* ap, double* x, int incx)
{
	//
	//     Test the input parameters.
	//
	int info = 0;
	if (!cblas_lsame(uplo, 'U') && !cblas_lsame(uplo, 'L')) {
		info = 1;
	}
	else if (!cblas_lsame(trans, 'N') &&
		!cblas_lsame(trans, 'T') &&
		!cblas_lsame(trans, 'C')) {
		info = 2;
	}
	else if (!cblas_lsame(diag, 'U') &&
		!cblas_lsame(diag, 'N')) {
		info = 3;
	}
	else if (n < 0) {
		info = 4;
	}
	else if (incx == 0) {
		info = 7;
	}
	if (info != 0) {
		cblas_xerbla("DTPSV ", info);
		return;
	}
	//
	//     Quick return if possible.
	//
	if (n == 0) {
		return;
	}

	//
	int i, ix, j, jx, k, kk, kx;
	double temp;
	double zero = 0.0;

	bool nounit = cblas_lsame(diag, 'N');
	//
	//     Set up the start point in X if the increment is not unity. This
	//     will be  ( N - 1 )*INCX  too small for descending loops.
	//
	if (incx <= 0) {
		kx = 1 - (n - 1) * incx;
	}
	else if (incx != 1) {
		kx = 0;
	}


	//
	//     Start the operations. In this version the elements of AP are
	//     accessed sequentially with one pass through AP.
	//
	if (cblas_lsame(trans, 'N')) {
		//
		//        Form  x := inv( A )*x.
		//
		if (cblas_lsame(uplo, 'U')) {
			//
			//  [ a0, a1, a2, a3, ]
			//  [     a4, a5, a6, ]
			//  [         a7, a8, ]
			//  [             a9  ]  
			//  
			//  x[3] = x[3] / a9  
			//  x[2] = x[2] - a8 * x[3]
			//  x[1] = x[1] - a6 * x[3] 
			//  x[0] = x[0] - a3 * x[3]
			// 
			//  x[2] = x[2] / a7 
			//  x[1] = x[1] - a5 * x[2] 
			//  x[0] = x[0] - a2 * x[2]
			// 
			//  x[1] = x[1] / a4 
			// 
			kk = (n * (n + 1)) / 2 - 1;
			if (incx == 1) {
				for (j = n - 1; j >= 0; j--) {
					if (x[j] != zero) {
						if (nounit) {
							x[j] = x[j] / ap[kk];
						}
						temp = x[j];
						k = kk;
						for (i = j - 1; i >= 0; i--)
						{
							k = k - (n - i - 1);

							//
							// printf("  x[%i] = x[%i] - temp * ap[%i]  \n", i, i, k);
							//

							x[i] = x[i] - temp * ap[k];
							// k = k - 1;  
						}
					}
					//kk = kk - j;
					kk = kk - (n - j + 1);
				}
			}
			else {
				jx = kx + (n - 1) * incx;
				for (j = n - 1; j >= 0; j--) {
					if (x[jx] != zero) {
						if (nounit) {
							x[jx] = x[jx] / ap[kk];
						}
						temp = x[jx];
						ix = jx;
						k = kk;
						//for (k = kk - 1; k < kk - j + 1; k--) 
						for (i = j - 1; i >= 0; i--)
						{
							k = k - (n - i - 1);
							ix = ix - incx;
							x[ix] = x[ix] - temp * ap[k];
						}
					}
					jx = jx - incx;
					kk = kk - (n - j + 1);
				}
			}
		}
		else {
			// uplo -- 'L'
			//
			//  [ a0,            ]
			//  [ a1, a2,        ]
			//  [ a3, a4, a5     ]
			//  [ a6, a7, a8, a9 ]  
			//
			//  
			//  x[0] = x[0] / a0  
			//  x[1] = x[1] - a1 * x[0]
			//  x[2] = x[2] - a3 * x[0] 
			//  x[3] = x[3] - a6 * x[0]
			// 
			// 
			kk = 0;
			if (incx == 1) {
				for (j = 0; j < n; j++) {
					if (x[j] != zero) {
						if (nounit) {
							x[j] = x[j] / ap[kk];
						}
						temp = x[j];
						k = kk;
						for (i = j + 1; i <= n - 1; i++) {
							k = k + i;
							//
							// printf("  x[%i] = x[%i] - temp * ap[%i]", i, i, k);
							//

							x[i] = x[i] - temp * ap[k];

							//
							// printf( " = %5.3f \n", x[i]);
							//
						}
					}
					kk = kk + (j + 2);
				}
			}
			else {
				jx = kx;
				for (j = 0; j < n; j++) {
					if (x[jx] != zero) {
						if (nounit) {
							x[jx] = x[jx] / ap[kk];
						}
						temp = x[jx];
						ix = jx;
						k = kk;
						//for (k = kk + 1; k < kk + n - j; k++) 
						for (i = j + 1; i <= n - 1; i++)
						{
							k = k + i;
							ix = ix + incx;
							x[ix] = x[ix] - temp * ap[k];
						}
					}
					jx = jx + incx;
					//kk = kk + (n - j + 1);
					kk = kk + (j + 2);
				}
			}
		}
	}
	else {
		//
		//        Form  x := inv( A**T )*x.
		//
		if (cblas_lsame(uplo, 'U')) {
			//
			//  Upper Triangle 
			// 
			//  [ a0, a1, a2, a3, ]
			//  [     a4, a5, a6, ]
			//  [         a7, a8, ]
			//  [             a9  ]  
			//  
			//  After transpose, the matrix is
			// 
			//  [ a0,            ]
			//  [ a1, a4,        ]
			//  [ a2, a5, a7,    ]
			//  [ a3, a6, a8, a9 ]  
			//   
			//
			//  
			//  x[0] = x[0] / a0  
			//  x[1] = x[1] - a1 * x[0]
			//  x[2] = x[2] - a2 * x[0] 
			//  x[3] = x[3] - a3 * x[0]
			// 
			// 


			kk = 0;
			if (incx == 1) {
				for (j = 0; j < n; j++) {
					temp = x[j];
					if (nounit) {
						temp = temp / ap[kk];
					}
					x[j] = temp;

					k = kk + 1;
					for (i = j + 1; i < n; i++) {
						x[i] = x[i] - temp * ap[k];
						//temp = temp - ap[k] * x[i];
						k = k + 1;
					}


					kk = kk + (n - j);
				}
			}
			else {
				jx = kx;
				for (j = 0; j < n; j++) {
					temp = x[jx];
					if (nounit) {
						temp = temp / ap[kk];
					}

					ix = kx + incx;
					k = kk + 1;
					for (i = j + 1; i < n; i++) {
						x[ix] = x[ix] - temp * ap[k];
						//temp = temp - ap[k] * x[i];
						ix = ix + incx;
						k = k + 1;
					}
					//for (k = kk; k < kk + j - 2; k++) {
					//	temp = temp - ap[k] * x[ix];
					//	ix = ix + incx;
					//}

					x[jx] = temp;
					jx = jx + incx;
					kk = kk + (n - j);

					kx = kx + incx;
				}
			}
		}
		else {
			//
			//  Lower Triangle 
			// 
			//  [ a0,         ]
			//  [ a1, a2,     ]
			//  [ a3, a4, a5, ]
			//  [ a6, ...     ]  
			//  
			//  After transpose, the matrix is
			//  [ a0, a1, a3, a6, ]    [x0] 
			//  [     a2, a4, a7, ]    [x1] 
			//  [         a5, a8, ]    [x2] 
			//  [             a9  ]    [x3]   
			// 
			//  
			//  x[3] = x[3] / a9  
			//  x[2] = x[2] - a8 * x[3]
			//  x[1] = x[1] - a7 * x[3] 
			//  x[0] = x[0] - a6 * x[3]
			// 
			//  x[2] = x[2] / a5 
			//  x[1] = x[1] - a4 * x[2] 
			//  x[0] = x[0] - a3 * x[2]
			// 
			//  x[1] = x[1] / a2 
			// 
			kk = (n * (n + 1)) / 2 - 1;
			if (incx == 1) {
				for (j = n - 1; j >= 0; j--) {
					temp = x[j];
					if (nounit) {
						temp = temp / ap[kk];
					}
					x[j] = temp;
					k = kk;
					for (i = j - 1; i >= 0; i--) {
						k = k - 1;
						x[i] = x[i] - temp * ap[k];
						//temp = temp - ap[k] * x[i]; 
					}
					kk = kk - (j + 1);
				}
			}
			else {
				kx = kx + (n - 1) * incx;
				jx = kx;
				for (j = n - 1; j >= 0; j--) {
					temp = x[jx];
					if (nounit) {
						temp = temp / ap[kk];
					}
					x[jx] = temp;
					ix = kx;
					k = kk;
					for (i = j - 1; i >= 0; i--) {
						k = k - 1;
						ix = ix - incx;
						x[ix] = x[ix] - temp * ap[k];
						//temp = temp - ap[k] * x[i]; 
					}
					kk = kk - (j + 1);
					jx = jx - incx;
					kx = kx - incx;
					//for (k = kk; k >= kk - (n - (j + 1)); k--) {
					//	temp = temp - ap[k] * x[ix];
					//	ix = ix - incx;
					//}
					//if (nounit) { 
					//	temp = temp / ap[kk - n + j]; 
					//} 
					//jx = jx - incx;
					//kk = kk - (n - j + 1);
				}
			}
		}
	}
	//
	return;
}
