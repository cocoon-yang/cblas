/*
 * cblas.cc
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


double cblas_dasum(const int N, const double *X, const int incX)
{
	double dtemp = 0.0;
	if ((N <= 0) || (incX <= 0))
	{
		printf("Error: in cblas_dasum(), Array number %6d less than ZERO.", N);
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
	//   code for m equal to 0
	//----------------------------------
	for ( i = m; i < N; i += 6)
	{
		dtemp += fabs(X[i]) + fabs(X[i + 1]) + fabs(X[i + 2]) + fabs(
				X[i + 3]) + fabs(X[i + 4]) + fabs(X[i + 5]);
	}
	return dtemp;
}



/*
*
* cblas_drot.c
*
* The program is a C version of drot.
*
*  Purpose
*  =======
*
*     applies a plane rotation.
*
*
*     Written by Chunfeng Yang
*     Date:  05/11/2014
*
*
*     Test:  Chunfeng Yang
*     Date: 05 November, 2014
*
*
*/
void cblas_drot(const int N, double *DX, const int INCX, double *DY,
	const int INCY, const double C, const double S)
{
	if (N <= 0)
	{
		return;
	}

	double DTEMP;
	int IX;
	int IY;
	int I;
	if ((INCX != 1) || (INCY != 1))
	{
		/*
		*       code for unequal increments or equal increments not equal
		*         to 1
		*/
		IX = 1;
		IY = 1;
		if (INCX < 0)
		{
			IX = (-N + 1) * INCX + 1;
		}
		if (INCY < 0)
		{
			IY = (-N + 1) * INCY + 1;
		}
		for (I = 0; I < N; I++)
		{
			DTEMP = C * DX[IX] + S * DY[IY];
			DY[IY] = C * DY[IY] - S * DX[I];
			DX[IX] = DTEMP;
			IX = IX + INCX;
			IY = IY + INCY;
		}
		return;
	}
	/*
	*       code for both increments equal to 1
	*/
	for (I = 0; I < N; I++)
	{
		DTEMP = C * DX[I] + S * DY[I];
		DY[I] = C * DY[I] - S * DX[I];
		DX[I] = DTEMP;
	}
	return;
}

void cblas_drotg(double* DA, double* DB, double C, double S)
{
	//    DOUBLE PRECISION R,ROE,SCALE,Z
	double R, ROE, SCALE, Z;


	//    ROE = DB
	ROE = *DB;

	//    IF (DABS(DA).GT.DABS(DB)) ROE = DA
	//    SCALE = DABS(DA) + DABS(DB)

	if (fabs(*DA) > fabs(*DB))
	{
		ROE = *DA;
	}

	//    IF (SCALE.NE.0.0d0) GO TO 10

	if (SCALE != 0)
	{
		goto Lab_10;
	}
	//    C = 1.0d0
	//    S = 0.0d0
	//    R = 0.0d0
	//    Z = 0.0d0

	C = 1.0;
	S = 0.0;
	R = 0.0;
	Z = 0.0;


	//    GO TO 20
	goto Lab_20;

	// 10 R = SCALE*DSQRT((DA/SCALE)**2+ (DB/SCALE)**2)
	//    R = DSIGN(1.0d0,ROE)*R
	//    C = DA/R
	//    S = DB/R
	//    Z = 1.0d0
	//    IF (DABS(DA).GT.DABS(DB)) Z = S
	//    IF (DABS(DB) >= DABS(DA) .AND. C.NE.0.0d0) Z = 1.0d0/C

Lab_10:
	R = SCALE*sqrt(pow((*DA / SCALE), 2) + pow((*DB / SCALE), 2));
	//R = DSIGN(1.0d0,ROE)*R;
	if (ROE > 0)
	{
		R = fabs(R);
	}
	else
	{
		R = -fabs(R);
	}

	C = *DA / R;
	S = *DB / R;
	Z = 1.0;
	if (fabs(*DA)>fabs(*DB))
	{
		Z = S;
	}
	if ((fabs(*DB) >= fabs(*DA)) && (C != 0.0))
	{
		Z = 1.0 / C;
	}

	// 20 DA = R
	//    DB = Z
	//    RETURN
	//    END

Lab_20:
	*DA = R;
	*DB = Z;

	return;
}


/*
*
* cblas_drotm.c
*
* The program is a C version of drot.
*
*  Purpose
*  =======
*
*     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
*
*     (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
*     (DY**T)
*
*     DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX  >=  0, ELSE
*     LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
*     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
*
*     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
*
*       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
*     H=(          )    (          )    (          )    (          )
*       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
*     SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         number of elements in input vector(s)
*
*  DX     (input/output) DOUBLE PRECISION array, dimension N
*         double precision vector with 5 elements
*
*  INCX   (input) INTEGER
*         storage spacing between elements of DX
*
*  DY     (input/output) DOUBLE PRECISION array, dimension N
*         double precision vector with N elements
*
*  INCY   (input) INTEGER
*         storage spacing between elements of DY
*
*  DPARAM (input/output)  DOUBLE PRECISION array, dimension 5
*     DPARAM(1)=DFLAG
*     DPARAM(2)=DH11
*     DPARAM(3)=DH21
*     DPARAM(4)=DH12
*     DPARAM(5)=DH22
*
*
*
*     Written by Chunfeng Yang
*     Date:  05/11/2014
*
*
*     Test:  Chunfeng Yang
*     Date: 05 November, 2014
*
*
*/

//      SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM)
//*     .. Scalar Arguments ..
//      INTEGER INCX,INCY,N
//*     ..
//*     .. Array Arguments ..
//      DOUBLE PRECISION DPARAM(5),DX(1),DY(1)

void cblas_drotm(const int N, double *DX, const int INCX, double *DY,
	const int INCY, const double* DPARAM)
{

	//	 *     .. Local Scalars ..
	//	       DOUBLE PRECISION DFLAG,DH11,DH12,DH21,DH22,TWO,W,Z,ZERO
	//	       INTEGER I,KX,KY,NSTEPS

	double DFLAG, DH11, DH12, DH21, DH22, W, Z;
	int I, KX, KY, NSTEPS;

	//    *     .. Data statements ..
	//          DATA ZERO,TWO/0.D0,2.D0/
	double ZERO = 0.0;
	double TWO = 2.0;

	//    DFLAG = DPARAM(1)
	//    IF (N.LE.0 .OR. (DFLAG+TWO.EQ.ZERO)) GO TO 140
	//    IF (.NOT. (INCX.EQ.INCY.AND.INCX.GT.0)) GO TO 70

	DFLAG = DPARAM[1];
	if (N <= 0 || DFLAG + TWO == ZERO)
	{
		goto L140;
	}
	if (!(INCX == INCY && INCX > 0))
	{
		goto L70;
	}

	NSTEPS = N * INCX;

	//       IF (DFLAG) 50,10,30


	if (DFLAG < 0.)
	{
		goto L50;
	}
	else if (DFLAG == 0)
	{
		goto L10;
	}
	else
	{
		goto L30;
	}

	//10 CONTINUE
	//DH12 = DPARAM(4)
	//DH21 = DPARAM(3)
	//DO 20 I = 1,NSTEPS,INCX
	//    W = DX(I)
	//    Z = DY(I)
	//    DX(I) = W + Z*DH12
	//    DY(I) = W*DH21 + Z
	//20 CONTINUE
	//    GO TO 140
L10: DH12 = DPARAM[4];
	DH21 = DPARAM[3];

	for (I = 1; I < NSTEPS; I += INCX)
	{
		W = DX[I];
		Z = DY[I];
		DX[I] = W + Z * DH12;
		DY[I] = W * DH21 + Z;
	}
	goto L140;

	//    30 CONTINUE
	//       DH11 = DPARAM(2)
	//       DH22 = DPARAM(5)
	//       DO 40 I = 1,NSTEPS,INCX
	//           W = DX(I)
	//           Z = DY(I)
	//           DX(I) = W*DH11 + Z
	//           DY(I) = -W + DH22*Z
	//    40 CONTINUE
	//       GO TO 140

L30: DH11 = DPARAM[2];
	DH22 = DPARAM[5];
	for (I = 1; I <= NSTEPS; I += INCX)
	{
		W = DX[I];
		Z = DY[I];
		DX[I] = W * DH11 + Z;
		DY[I] = -W + DH22 * Z;
	}
	goto L140;

	//	   50 CONTINUE
	//	      DH11 = DPARAM(2)
	//	      DH12 = DPARAM(4)
	//	      DH21 = DPARAM(3)
	//	      DH22 = DPARAM(5)
	//	      DO 60 I = 1,NSTEPS,INCX
	//	          W = DX(I)
	//	          Z = DY(I)
	//	          DX(I) = W*DH11 + Z*DH12
	//	          DY(I) = W*DH21 + Z*DH22
	//	   60 CONTINUE
	//	      GO TO 140

L50:
	DH11 = DPARAM[2];
	DH12 = DPARAM[4];
	DH21 = DPARAM[3];
	DH22 = DPARAM[5];
	for (I = 1; I <= NSTEPS; I += INCX)
	{
		W = DX[I];
		Z = DY[I];
		DX[I] = W * DH11 + Z * DH12;
		DY[I] = W * DH21 + Z * DH22;
	}
	goto L140;


	//	   70 CONTINUE
	//	      KX = 1
	//	      KY = 1
	//	      IF (INCX.LT.0) KX = 1 + (1-N)*INCX
	//	      IF (INCY.LT.0) KY = 1 + (1-N)*INCY

L70:
	KX = 1;
	KY = 1;
	if (INCX < 0)
	{
		KX = 1 + (1 - N) * INCX;
	}
	if (INCY < 0)
	{
		KY = 1 + (1 - N) * INCY;
	}


	//    IF (DFLAG) 120,80,100
	if (DFLAG < 0.)
	{
		goto L120;
	}
	else if (DFLAG == 0)
	{
		goto L80;
	}
	else
	{
		goto L100;
	}


	//	   80 CONTINUE
	//	      DH12 = DPARAM(4)
	//	      DH21 = DPARAM(3)
	//	      DO 90 I = 1,N
	//	          W = DX(KX)
	//	          Z = DY(KY)
	//	          DX(KX) = W + Z*DH12
	//	          DY(KY) = W*DH21 + Z
	//	          KX = KX + INCX
	//	          KY = KY + INCY
	//	   90 CONTINUE
	//	      GO TO 140
L80:
	DH12 = DPARAM[4];
	DH21 = DPARAM[3];
	for (I = 1; I <= N; I++)
	{
		W = DX[KX];
		Z = DY[KY];
		DX[KX] = W + Z * DH12;
		DY[KY] = W * DH21 + Z;
		KX = KX + INCX;
		KY = KY + INCY;
	}
	goto L140;



	//	  100 CONTINUE
	//	      DH11 = DPARAM(2)
	//	      DH22 = DPARAM(5)
	//	      DO 110 I = 1,N
	//	          W = DX(KX)
	//	          Z = DY(KY)
	//	          DX(KX) = W*DH11 + Z
	//	          DY(KY) = -W + DH22*Z
	//	          KX = KX + INCX
	//	          KY = KY + INCY
	//	  110 CONTINUE
	//	      GO TO 140

L100:
	DH11 = DPARAM[2];
	DH22 = DPARAM[5];
	for (I = 1; I <= N; I++)
	{
		W = DX[KX];
		Z = DY[KY];
		DX[KX] = W*DH11 + Z;
		DY[KY] = -W + DH22*Z;
		KX = KX + INCX;
		KY = KY + INCY;
	}
	goto L140;


	//	  120 CONTINUE
	//	      DH11 = DPARAM(2)
	//	      DH12 = DPARAM(4)
	//	      DH21 = DPARAM(3)
	//	      DH22 = DPARAM(5)
	//	      DO 130 I = 1,N
	//	          W = DX(KX)
	//	          Z = DY(KY)
	//	          DX(KX) = W*DH11 + Z*DH12
	//	          DY(KY) = W*DH21 + Z*DH22
	//	          KX = KX + INCX
	//	          KY = KY + INCY
	//	  130 CONTINUE

L120:
	DH11 = DPARAM[2];
	DH12 = DPARAM[4];
	DH21 = DPARAM[3];
	DH22 = DPARAM[5];
	for (I = 1; I <= N; I++)
	{
		W = DX[KX];
		Z = DY[KY];
		DX[KX] = W * DH11 + Z * DH12;
		DY[KY] = W * DH21 + Z * DH22;
		KX = KX + INCX;
		KY = KY + INCY;
	}
	goto L140;

L140:
	return;
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





/*  Purpose
*  =======
**
*     scales a vector by a constant.
*     uses unrolled loops for increment equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*/
void cblas_dscal(const int N, const double DA, double *DX, const int INCX)
{

	if ((N <= 0) || (INCX <= 0))
	{
		return;
	}

	if (INCX != 1)
	{
		/*
		*        code for increment not equal to 1
		*/
		int NINCX = N * INCX;
		for (int I = 0; I < NINCX; I += INCX)
		{
			DX[I] = DA * DX[I];

		}
		return;
	}
	/*
	*        code for increment equal to 1
	*
	*
	*        clean-up loop
	*/
	int M = (N % 5);
	if (M != 0)
	{
		for (int I = 0; I < M; I++)
		{
			DX[I] = DA * DX[I];
		}

		if (N < 5)
		{
			return;
		}
	}
	// MP1 = M + 1
	int MP1 = M;
	for (int I = MP1; I < N; I += 5)
	{
		DX[I] = DA * DX[I];
		DX[I + 1] = DA * DX[I + 1];
		DX[I + 2] = DA * DX[I + 2];
		DX[I + 3] = DA * DX[I + 3];
		DX[I + 4] = DA * DX[I + 4];
	}

	return;
}

void cblas_sscal(const int N, const float ALPHA, float *X, const int INCX)
{

	if ((N <= 0) || (INCX <= 0))
	{
		return;
	}

	if (INCX != 1)
	{
		/*
		*        code for increment not equal to 1
		*/
		int NINCX = N * INCX;
		for (int I = 0; I < NINCX; I += INCX)
		{
			X[I] = ALPHA * X[I];

		}
		return;
	}
	/*
	*        code for increment equal to 1
	*
	*
	*        clean-up loop
	*/
	int M = (N % 5);
	if (M != 0)
	{
		for (int I = 0; I < M; I++)
		{
			X[I] = ALPHA * X[I];
		}

		if (N < 5)
		{
			return;
		}
	}
	// MP1 = M + 1
	int MP1 = M;
	for (int I = MP1; I < N; I += 5)
	{
		X[I] = ALPHA * X[I];
		X[I + 1] = ALPHA * X[I + 1];
		X[I + 2] = ALPHA * X[I + 2];
		X[I + 3] = ALPHA * X[I + 3];
		X[I + 4] = ALPHA * X[I + 4];
	}

	return;
}

/*
*
*
*  Purpose
*  =======
*
*     interchanges two vectors.
*         X <--> Y
*     uses unrolled loops for increments equal one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*/
//      SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
void cblas_dswap(const int N, double *DX, const int INCX, double *DY,
	const int INCY)
{

	//     .. Local Scalars ..
	double DTEMP;
	int I, IX, IY, M, MP1;
	/*     ..
	*     .. Intrinsic Functions ..
	*/
	// INTRINSIC MOD;//  ???

	if (N <= 0)
	{
		return;
	}
	//if(INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
	if ((INCX != 1) || (INCY != 1))
	{
		/*
		*       code for unequal increments or equal increments not equal
		*         to 1
		*/
		IX = 1;
		IY = 1;
		if (INCX < 0)
			IX = (-N + 1) * INCX + 1;
		if (INCY < 0)
			IY = (-N + 1) * INCY + 1;
		for (I = 1; I < N; I++)
		{
			DTEMP = DX[IX];
			DX[IX] = DY[IY];
			DY[IY] = DTEMP;
			IX = IX + INCX;
			IY = IY + INCY;
		}
		return;
	}
	/*
	*       code for both increments equal to 1
	*
	*
	*       clean-up loop
	*/
	M = N % 3;
	if (M != 0)
	{
		for (I = 1; I < M; I++)
		{
			DTEMP = DX[I];
			DX[I] = DY[I];
			DY[I] = DTEMP;
		}
		if (N < 3)
			return;
	}
	// MP1 = M + 1
	MP1 = M;
	for (I = MP1; I < N; I += 3)
	{
		DTEMP = DX[I];
		DX[I] = DY[I];
		DY[I] = DTEMP;
		DTEMP = DX[I + 1];
		DX[I + 1] = DY[I + 1];
		DY[I + 1] = DTEMP;
		DTEMP = DX[I + 2];
		DX[I + 2] = DY[I + 2];
		DY[I + 2] = DTEMP;
	}
	return;
}


/*
*  Purpose
*  =======
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*
*  -- This version written on 25-October-1982.
*     Modified on 14-October-1993 to inline the call to DLASSQ.
*     Sven Hammarling, Nag Ltd.
*
*/
double cblas_dnrm2(const int N, const double *X, const int incX)
{
	double ONE, ZERO;
	ONE = 1.0;
	ZERO = 0.0;

	//    .. Local Scalars ..
	double ABSXI, NORM, SCALE, SSQ;
	int IX;

	if (N < 1 || incX < 1)
	{
		NORM = ZERO;
	}
	else if (N == 1)
	{
		NORM = fabs(X[0]);
	}
	else
	{
		SCALE = ONE;
		SSQ = ZERO;
		/*        The following loop is equivalent to this call to the LAPACK
		*        auxiliary routine:
		*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
		*/
		for (IX = 0; IX < (1 + (N - 1) * incX); IX += incX)
		{
			if (X[IX] != ZERO)
			{
				ABSXI = fabs(X[IX]);

				if (SCALE>ABSXI)
				{
					SSQ = ONE + SSQ * (SCALE / ABSXI) * (SCALE / ABSXI);
					SCALE = ABSXI;
				}
				else
				{
					SSQ = SSQ + (ABSXI / SCALE) * (ABSXI / SCALE);
				}
			}
		}
		NORM = SCALE * sqrt(SSQ);
	} 
	return NORM;
}



/***
DLASWP performs a series of row interchanges on the matrix A.
One row interchange is initiated for each of rows K1 through K2 of A.

M: int, The number of rows of the matrix A.  M >= 0.
N: int, The number of columns of the matrix A.  N >= 0.
A: double*,  A is DOUBLE PRECISION array, dimension (LDA,N)
On entry, the matrix of column dimension N to which the row
interchanges will be applied.
On exit, the permuted matrix.
LDA: int, The leading dimension of the array A.  LDA >= max(1,M).
K1: int, The first element of IPIV for which a row interchange will
be done.
K2: int, (K2-K1+1) is the number of elements of IPIV for which a row
interchange will be done.
IPIV: int*,  IPIV is INTEGER array, dimension (K1+(K2-K1)*abs(INCX))
The vector of pivot indices. Only the elements in positions
K1 through K1+(K2-K1)*abs(INCX) of IPIV are accessed.
IPIV(K1+(K-K1)*abs(INCX)) = L implies rows K and L are to be
interchanged.
INCX: int, The increment between successive values of IPIV. If INCX
is negative, the pivots are applied in reverse order.
***/
void cblas_dlaswp(int n, double* a, int lda, int k1, int k2, int* ipiv, int incx)
{
	int  i, i1, i2, inc, ip, ix, ix0, j, k, n32;
	double   temp;

	//      ..
	//      .. Executable Statements ..
	// 
	//      Interchange row I with row IPIV(K1+(I-K1)*abs(INCX)) for each of rows
	//      K1 through K2.
	// 
	if (incx > 0) {
		ix0 = k1;
		i1 = k1;
		i2 = k2;
		inc = 1;
	}
	else if (incx < 0) {
		ix0 = k1 + (k1 - k2)*incx;
		i1 = k2;
		i2 = k1;
		inc = -1;
	}
	else {
		return;
	}
	// 
	n32 = (n / 32) * 32;
	if (n32 != 0) {
		for (j = 0; j < n32; j += 32) {
			ix = ix0;
			for (i = i1; i < i2; i += inc) {
				ip = ipiv[ix];
				if (ip != i) {
					for (k = j; k < j + 31; k++) {
						temp = a[i*lda + k];
						a[i*lda + k] = a[ip*lda + k];
						a[ip*lda + k] = temp;
					}
				}
				ix = ix + incx;
			}
		}
	}
	if (n32 != n) {
		n32 = n32 + 1;
		ix = ix0;
		for (i = i1; i < i2; i += inc) {
			ip = ipiv[ix];
			if (ip != i) {
				for (k = n32; k < n; k++) {
					temp = a[i*lda + k];
					a[i*lda + k] = a[ip*lda + k];
					a[ip*lda + k] = temp;
				}
			}
			ix = ix + incx;
		}
	}
	// 
	return;


}