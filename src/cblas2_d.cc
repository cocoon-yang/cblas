#include "cblas.h"
#include <utility>
#include "math.h" 


//==============================================================
// Level 2

/*
Purpose:

SGER computes A := alpha*x*y' + A.

Discussion:

SGER performs the rank 1 operation

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
				TEMP = ALPHA * X[I];
				for (J = 0; J < N; J = J + INCY)
				{
					A[(I)* LDA + (J)] = A[(I)* LDA + (J)] + Y[J] * TEMP;
				}
			}
			//JY = JY + INCY;
		}
	}
	else
	{
		for (I = 0; I < M; I += INCX)
		{
			if (X[I] != ZERO)
			{
				TEMP = ALPHA * X[I];
				for (J = 0; J < N; J = J + INCY)
				{
					A[(I)* LDA + (J)] = A[(I)* LDA + (J)] + Y[J] * TEMP;
				}
			}
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
int cblas_dgemv(char TRANS, int M, int N, double ALPHA, double* A, int LDA,
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
	else if (LDA < std::max(1, M))
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
	//	      IF (LSAME(TRANS,'N')) THEN
	if (cblas_lsame(TRANS, 'N'))
	{
		/*
		*        Form  y := alpha*A*x + y.
		*/
		// printf("Form  y := alpha*A*x + y.  \n");
		JY = KY;
		if (INCX == 1)
		{
			for (J = 0; J < N; J++) // C/C++ version
									//for (J = 1; J < N; J++)  // Fortran version
			{
				TEMP = ZERO;
				for (I = 0; I < M; I++)
					// for (I = 1; I < M; I++) // Fortran version
				{
					TEMP = TEMP + A[(I)+(J)* LDA] * X[I];
				}
				Y[JY] = Y[JY] + ALPHA * TEMP;
				JY = JY + INCY;
			}
		}
		else
		{
			for (J = 0; J < N; J++) // C/C++ version
									//for (J = 1; J < N; J++)  // Fortran version
			{
				TEMP = ZERO;
				IX = KX;
				for (I = 0; I < M; I++)
					// for (I = 1; I < M; I++) // Fortran version
				{
					TEMP = TEMP + A[(I)+(J)* LDA] * X[IX];
					IX = IX + INCX;
				}
				Y[JY] = Y[JY] + ALPHA * TEMP;
				JY = JY + INCY;
			}
		}
	}
	else
	{
		//	     *
		//	     *        Form  y := alpha*A'*x + y.
		//	     *
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
						Y[I] = Y[I] + TEMP * A[I + J * LDA];
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
						Y[IY] = Y[IY] + TEMP * A[(I)+(J)* LDA];
						IY = IY + INCY;
					}
				}
				JX = JX + INCX;
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
void dtrmv(char uplo, char trans, char diag, int n, double a[], int lda,
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
		xerbla("DTRMV", info);
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
				for (j = n - 1; 0 <= j; j--)
				{
					temp = x[j];
					if (nounit)
					{
						temp = temp * a[j + j*lda];
					}
					for (i = j - 1; 0 <= i; i--)
					{
						temp = temp + a[i + j*lda] * x[i];
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
						temp = temp * a[j + j*lda];
					}
					for (i = j - 1; 0 <= i; i--)
					{
						ix = ix - incx;
						temp = temp + a[i + j*lda] * x[ix];
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
						temp = temp * a[j + j*lda];
					}
					for (i = j + 1; i < n; i++)
					{
						temp = temp + a[i + j*lda] * x[i];
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
						temp = temp * a[j + j*lda];
					}
					for (i = j + 1; i < n; i++)
					{
						ix = ix + incx;
						temp = temp + a[i + j*lda] * x[ix];
					}
					x[jx] = temp;
					jx = jx + incx;
				}
			}
		}
	}
	//
	//  Form x := A'*x.
	//
	else
	{
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
							x[i] = x[i] + temp * a[i + j*lda];
						}
						if (nounit)
						{
							x[j] = x[j] * a[j + j*lda];
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
							x[ix] = x[ix] + temp * a[i + j*lda];
							ix = ix + incx;
						}
						if (nounit)
						{
							x[jx] = x[jx] * a[j + j*lda];
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
							x[i] = x[i] + temp * a[i + j*lda];
						}
						if (nounit)
						{
							x[j] = x[j] * a[j + j*lda];
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
							x[ix] = x[ix] + temp * a[i + j*lda];
							ix = ix - incx;
						}
						if (nounit)
						{
							x[jx] = x[jx] * a[j + j*lda];
						}
					}
					jx = jx - incx;
				}
			}
		}
	}

	return;
}



