/*
 * clas.cc
 *
 *  Created on: 2018-12-24
 *      Author: Yang
 */
#include "cblas.h"
#include <utility>
#include "math.h"

//==============================================================
// Level 1

/*
 * cblas_dasum.c
 *
 * The program is a C version of dasum.
 * It calls the fortran wrapper before calling dasum.
 *
 * Written by Yang Chunfeng.  27/05/2012
 *
 *
 *  Test
 *    by Yang
 *    Date: 09 Feb, 2013
 *
 *
 *  COPYRIGHT STATEMENT
 *
 *  You may only use this program for your own private purposes.
 *  You are not allowed, in any circumstances, to distribute this
 *  program (including its source code, executable and any other
 *  files related to it, either in their original version or any
 *  modifications introduced by you, the authors or any other party)
 *  in whole or in part, either freely or otherwise, in any medium,
 *  without the prior written consent of the copyright holders.
 *
 *  DISCLAIMER
 *
 *  This program (including its source code, executable and any
 *  other files related to it) is provided "as is" without warranty
 *  of any kind, either expressed or implied, including, but not
 *  limited to, any implied warranties of fitness for purpose.
 *  In particular, THIS PROGRAM IS BY NO MEANS GUARANTEED TO BE FREE
 *  FROM ERRORS.
 *  The results produced by this program are in no way garanteed to
 *  be fit for any purpose.
 *  This program (or any modification incorporated to it by you, the
 *  authors or any other party) will be used entirely at your own
 *  risk.
 *  Under no circumstances will the authors/copyright holders be
 *  liable to anyone for damages, including any general, special,
 *  incidental or consequential damages arising from the use or
 *  inability to use the program (including, but not limited to,
 *  loss or corruption of data, failure of the program to operate in
 *  any particular way as well as damages arising from the use of
 *  any results produced by the program for any purpose).
 *
 *  CONDITIONS OF USE
 *
 *  You may only use this program if you fully understand and agree
 *  with the terms of the above disclaimer. You must not use this
 *  program if you do not agree with or do not understand (fully or
 *  in part) these conditions of use.
 *
 */



#include "math.h"
double cblas_dasum(const int N, const double *X, const int incX)
{
	double dtemp = 0.0;
	if ((N <= 0) || (incX <= 0))
	{
		printf("Error: in cblas_dasum(), Array number %6d less than zero.", N);
		return dtemp;
	}
	int i;
	if (incX != 1)
	{
		//--------------------------------------------
		//        code for increment not equal to 1
		//--------------------------------------------
		int nincX = N * incX;
		for ( i = 0; i < nincX; i += incX)
		{
			dtemp = dtemp + fabs(X[i]);
		}
		return dtemp;
	}
	//--------------------------------------------
	//        code for increment equal to 1
	//
	//
	//        clean-up loop
	//--------------------------------------------

	int m = N % 6;
	if (0 != m)
	{
		for ( i = 0; i < m; i++)
		{
			dtemp = dtemp + fabs(X[i]);
		}
		if (N < 6)
		{
			return dtemp;
		}
	}
	//----------------------------------
	//        code for m equal to 0
	//----------------------------------
	for ( i = m; i < N; i += 6)
	{
		dtemp += fabs(X[i]) + fabs(X[i + 1]) + fabs(X[i + 2]) + fabs(
				X[i + 3]) + fabs(X[i + 4]) + fabs(X[i + 5]);
	}
	return dtemp;
}


void cblas_daxpy(const int N, const double DA, const double *X,
		const int incX, double *Y, const int incY)
{


	if (N <= 0)
	{
		return;
	}
	if (DA == 0.0)
	{
		return;
	}
	//
	//        code for unequal increments or equal increments
	//          not equal to 1
	//
	int I;
	if ((incX != 1) || (incY != 1))
	{

		int ix = 1;
		int iy = 1;
		if (incX < 0)
		{
			ix = (1 - N) * incX + 1;
		}
		if (incY < 0)
		{
			iy = (1 - N) * incY + 1;
		}
		for ( I = 0; I < N; I++)
		{
			Y[iy] = Y[iy] + DA * X[ix];
			ix = ix + incX;
			iy = iy + incY;
		}
		return;
	}

	//        code for both increments equal to 1
	//        clean-up loop
	//double dx [N];
	//double dy [N];

	int m = N % 4; // mod(N, 4);
	if (0 != m)
	{
		for (I = 0; I < m; I++)
		{
			Y[I] = Y[I] + DA * X[I];
		}
		if (N < 4)
		{
			return;
		}
	}

	//40 mp1 = m + 1
	int mp1 = m;
	for (I = mp1; I < N; I += 4)
	{
		Y[I] = Y[I] + DA * X[I];
		Y[I + 1] = Y[I + 1] + DA * X[I + 1];
		Y[I + 2] = Y[I + 2] + DA * X[I + 2];
		Y[I + 3] = Y[I + 3] + DA * X[I + 3];
	}
	return;
}

void cblas_dcopy(const int N, const double *X, const int incX, double *Y,
		const int incY)
{

	if (N <= 0)
	{
		return;
	}
	//
	//        code for unequal increments or equal increments
	//          not equal to 1
	//
	int i;
	if ((incX != 1) || (incY != 1))
	{
		int ix = 1;
		int iy = 1;
		if (incX < 0)
		{
			ix = (1 - N) * incX + 1;
		}
		if (incY < 0)
		{
			iy = (1 - N) * incY + 1;
		}
		for ( i = 0; i < N; i++)
		{
			Y[iy] = X[ix];
			ix = ix + incX;
			iy = iy + incY;
		}
		return;
	}

	//
	//        code for both increments equal to 1
	//
	//
	//        clean-up loop
	//
	int m = (N % 7);
	if (0 != m)
	{
		for ( i = 0; i < m; i++)
		{
			Y[i] = X[i];
		}
		if (N < 7)
		{
			return;
		}
	}
	int mp1 = m ;
	for ( i = mp1; i < N; i += 7)
	{
		Y[i] = X[i];
		Y[i + 1] = X[i + 1];
		Y[i + 2] = X[i + 2];
		Y[i + 3] = X[i + 3];
		Y[i + 4] = X[i + 4];
		Y[i + 5] = X[i + 5];
		Y[i + 6] = X[i + 6];
	}
	return;

}


/*
 * cblas_ddot.c
 *
 * The program is a C version of ddot.
 */

double cblas_ddot(const int N, const double *X, const int incX,
		const double *Y, const int incY)
{
	double result = 0.0;
	if (N < 0)
	{
		return result;
	}

	//    IX = 1
	//    IY = 1
	//    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
	//    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
	//    DO 10 I = 1,N
	//        DTEMP = DTEMP + DX(IX)*DY(IY)
	//        IX = IX + INCX
	//        IY = IY + INCY
	// 10 CONTINUE
	//    DDOT = DTEMP
	//    RETURN
	//
	//        code for unequal increments or equal increments
	//          not equal to 1
	//
	int i;
	if ((incX != 1) || (incY != 1))
	{
		int ix = 1;
		int iy = 1;
		if (incX < 0)
		{
			ix = (1 - N) * incX + 1;
		}
		if (incY < 0)
		{
			iy = (1 - N) * incY + 1;
		}
		for ( i = 0; i < N; i++)
		{
			result += Y[iy] * X[ix];
			ix = ix + incX;
			iy = iy + incY;
		}
		return result;
	}
	//
	//        code for both increments equal to 1
	//
	//
	//        clean-up loop
	//

	int m = (N % 5);//int m = mod(N, 5);
	if (0 != m)
	{
		for ( i = 0; i < m; i++)
		{
			result += Y[i] * X[i];
		}
		if (N < 5)
		{
			return result;
		}
	}

	//int mp1 = m + 1;
	int mp1 = m ;
	for ( i = mp1; i < N; i += 5)
	{
		result += (Y[i] * X[i] + Y[i + 1] * X[i + 1] + Y[i + 2] * X[i + 2]
				+ Y[i + 3] * X[i + 3] + Y[i + 4] * X[i + 4]);

	}
	return result;

}





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
	else if (LDA < MAX(1, M))
	{
		INFO = 9;
	}
	if (INFO != 0)
	{
		cblas_xerbla("DGER ", INFO);
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
		//JY = 1; // Fortran version
	}
	else
	{
		JY = 0 - (N - 1) * INCY;
		//JY = 1 - (N - 1) * INCY; // Fortran version
	}

	if (INCX == 1)
	{
		for (J = 0; J < N; J++)
		//for (J = 1; J < N; J++)  // fortra
		{
			if (Y[JY] != ZERO)
			{
				TEMP = ALPHA * Y[JY];

				//for (I = 1; I < M; I++)  // fortra
				for (I = 0; I < M; I++)
				{
					A[(I) + (J) * LDA] = A[(I) + (J) * LDA] + X[I] * TEMP;
				}
			}
			JY = JY + INCY;
		}
	}
	else
	{
		if (INCX > 0)
		{
			//KX = 1;  // Fortran version
			KX = 0;
		}
		else
		{
			//KX = 1 - (M - 1) * INCX; // Fortran version
			KX = -(M - 1) * INCX;
		}
		for (J = 0; J < N; J++)
		//for (J = 1; J < N; J++) // Fortran version
		{
			if (Y[JY] != ZERO)
			{
				TEMP = ALPHA * Y[JY];
				IX = KX;
				for (I = 0; I < M; I++)
				//for (I = 1; I < M; I++)  // Fortran version
				{
					A[(I) + (J) * LDA] = A[(I) + (J) * LDA] + X[IX] * TEMP;
				}
				IX = IX + INCX;
			}
		}
		JY = JY + INCY;
	}
	//}

	return 0;

}


int cblas_dgemv(char TRANS, int M, int N, double ALPHA, double* A, int LDA,
		double* X, int INCX, double BETA, double* Y, int INCY)
{
	//     DOUBLE PRECISION ONE,ZERO
	//     PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
	//*     ..
	//*     .. Local Scalars ..
	//     DOUBLE PRECISION TEMP
	//     INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
	double ONE = 1.0;
	double ZERO = 0.0;

	double TEMP;
	int I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY;

	/*
	 *     Test the input parameters.
	 */
	//	      INFO = 0
	//	      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
	//	     +    .NOT.LSAME(TRANS,'C')) THEN
	//	          INFO = 1
	//	      ELSE IF (M.LT.0) THEN
	//	          INFO = 2
	//	      ELSE IF (N.LT.0) THEN
	//	          INFO = 3
	//	      ELSE IF (LDA.LT.MAX(1,M)) THEN
	//	          INFO = 6
	//	      ELSE IF (INCX.EQ.0) THEN
	//	          INFO = 8
	//	      ELSE IF (INCY.EQ.0) THEN
	//	          INFO = 11
	//	      END IF
	//	      IF (INFO.NE.0) THEN
	//	          CALL XERBLA('DGEMV ',INFO)
	//	          RETURN
	//	      END IF
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
	else if (LDA < MAX(1, M))
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
		cblas_xerbla("DGEMV ", INFO);
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
		JX = KX;
		if (INCY == 1)
		{
			for (J = 0; J < N; J++) // C/C++ version
			//for (J = 1; J < N; J++)  // Fortran version
			{
				if (X[JX] != ZERO)
				{
					TEMP = ALPHA * X[JX];
					for (I = 0; I < M; I++)
					// for (I = 1; I < M; I++) // Fortran version
					{
						Y[I] = Y[I] + TEMP * A[I + J * LDA];
					}
				}
				JX = JX + INCX;
			}
		}
		else
		{
			for (J = 0; J < N; J++) // C/C++ version
			//for (J = 1; J < N; J++)  // Fortran version
			{
				if (X[JX] != ZERO)
				{
					TEMP = ALPHA * X[JX];
					IY = KY;
					for (I = 0; I < M; I++)
					// for (I = 1; I < M; I++) // Fortran version
					{
						Y[IY] = Y[IY] + TEMP * A[(I) + (J) * LDA];
						IY = IY + INCY;
					}
				}
				JX = JX + INCX;
			}
		}
	}
	else
	{
		//	     *
		//	     *        Form  y := alpha*A'*x + y.
		//	     *
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
					TEMP = TEMP + A[(I) + (J) * LDA] * X[I];
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
					TEMP = TEMP + A[(I) + (J) * LDA] * X[IX];
					IX = IX + INCX;
				}
				Y[JY] = Y[JY] + ALPHA * TEMP;
				JY = JY + INCY;
			}
		}
	}

}


void xerbla( char* SRNAME, int INFO)
{
//	std::cout << "CBLAS: On entry to " << SRNAME << " parameter number " << INFO
//	<< " had an illegal value" << std::endl;
	printf( "CBLAS: On entry to  %s parameter number %c had an illegal value.", SRNAME, INFO);

}

void cblas_xerbla( char* SRNAME, int INFO)
{
//	std::cout << "CBLAS: On entry to " << SRNAME << " parameter number " << INFO
//	<< " had an illegal value" << std::endl;
	printf( "CBLAS: On entry to  %s parameter number %c had an illegal value.", SRNAME, INFO);

}


int MAX( int first, int second )
{
	return first? second: first > second;
}

bool cblas_lsame(char CA, char CB)
{
	/* System generated locals */
	bool ret_val;

	/* Local variables */
	//int inta, intb, zcode;
	int INTA, INTB, ZCODE;

	/*     Test if the characters are equal */

	ret_val = (CA == CB);
	if (ret_val)
	{
		return ret_val;
	}

	/*     Now test for equivalence if both characters are alphabetic. */

	ZCODE = 'Z';

	/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
	/*     machines, on which ICHAR returns a value with bit 8 set. */
	/*     ICHAR('A') on Prime machines returns 193 which is the same as */
	/*     ICHAR('A') on an EBCDIC machine. */

	INTA = (unsigned char) CA;
	INTB = (unsigned char) CB;

	if (ZCODE == 90 || ZCODE == 122)
	{

		/*        ASCII is assumed - ZCODE is the ASCII code of either lower or */
		/*        upper case 'Z'. */

		if (INTA >= 97 && INTA <= 122)
		{
			INTA += -32;
		}
		if (INTB >= 97 && INTB <= 122)
		{
			INTB += -32;
		}

	}
	else if (ZCODE == 233 || ZCODE == 169)
	{

		/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or */
		/*        upper case 'Z'. */

		if (INTA >= 129 && INTA <= 137 || INTA >= 145 && INTA <= 153 || INTA
				>= 162 && INTA <= 169)
		{
			INTA += 64;
		}
		if (INTB >= 129 && INTB <= 137 || INTB >= 145 && INTB <= 153 || INTB
				>= 162 && INTB <= 169)
		{
			INTB += 64;
		}

	}
	else if (ZCODE == 218 || ZCODE == 250)
	{

		/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code */
		/*        plus 128 of either lower or upper case 'Z'. */

		if (INTA >= 225 && INTA <= 250)
		{
			INTA += -32;
		}
		if (INTB >= 225 && INTB <= 250)
		{
			INTB += -32;
		}
	}
	ret_val = INTA == INTB;

	/*     RETURN */

	/*     End of LSAME */

	return ret_val;
}
