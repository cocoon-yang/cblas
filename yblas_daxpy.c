/*
 * cblas_daxpy.c
 *
 *    constant times a vector plus a vector.
 *    uses unrolled loops for increments equal to one.
 *
 * The program is a C version of daxpy.
 *
 * Written by Yang Chunfeng.
 *
 */
#include "yblas.h"
void cblas_daxpy(const int N, const double alpha, const double *X,
  	const int incX, double *Y, const int incY)
{
	if (N <= 0)
	{
		return;
	}
	if (alpha == 0.0)
	{
		return;
	}
	//
	//        code for unequal increments or equal increments
	//          not equal to 1
	//
	if ((incx != 1) || (incy != 1))
	{
		int ix = 1;
		int iy = 1;
		if (incX < 0)
		{
			ix = (1 - N) * incX + 1;
		}
		if (incy < 0)
		{
			iy = (1 - N) * incY + 1;
		}
		for (int i = 0; i < N; i++)
		{
			Y[iy] = Y[iy] + alpha * X[ix];
			ix = ix + incX;
			iy = iy + incY;
		}
		return;
	}

	//        code for both increments equal to 1
	//        clean-up loop

	int m = mod(N, 4);
	if (0 != m)
	{
		for (int i = 0; i < m; i++)
		{
			dy[i] = dy[i] + da * dx[i];
		}
		if (n < 4)
		{
			return;
		}
	}

	//40 mp1 = m + 1
	int mp1 = m;
	for (int i = mp1; i < n; i += 4)
	{
		dy[i] = dy[i] + da * dx[i];
		dy[i + 1] = dy[i + 1] + da * dx[i + 1];
		dy[i + 2] = dy[i + 2] + da * dx[i + 2];
		dy[i + 3] = dy[i + 3] + da * dx[i + 3];
	}
	return;
}



/*******************************************
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     constant times a vector plus a vector.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          DY(I) = DY(I) + DA*DX(I)
          DY(I+1) = DY(I+1) + DA*DX(I+1)
          DY(I+2) = DY(I+2) + DA*DX(I+2)
          DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
      END
****************************/
