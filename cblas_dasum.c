/*
 * yblas_dasum.c
 *
 *  Created on: Jul 6, 2013
 *      Author: yang
 */

/*
 *  Purpose
 *  =======
 *
 *     takes the sum of the absolute values.
 *
 */
void cblas_dasum(int N, double*DX, int INCX)
{

  //     .. Local Scalars ..
	double DTEMP;
	int I, M, MP1, NINCX;
	//     ..
	DASUM = 0.0;
	DTEMP = 0.0;
	if (N <= 0 || INCX <= 0)
	{
		return;
	}
	if (INCX != 1)
	{
		//
		//        code for increment not equal to 1
		//
		NINCX = N;//INCX
		for (int I = 0; i < NINCX; i += INCX)
		{
			DTEMP = DTEMP + fabs(DX[I]);
		}
		DASUM = DTEMP;
		return;
	}
	//
	//        code for increment equal to 1
	//
	//
	//        clean-up loop
	//
	M = MOD(N, 6);
	if (M != 0)
	{
		for (int I = 0; i < M; i++)
		{
			DTEMP = DTEMP + fabs(DX[I]);
		}
		if (N < 6)
		{
			DASUM = DTEMP;
			return;
		}
	}
	//MP1 = M + 1;
  	MP1 = M ;
	for (int I = MP1; i < N; i += 6)
	{
		DTEMP = DTEMP + fabs(DX[I]) + fabs(DX[I + 1]) + fabs(DX[I + 2]) + fabs(
				DX[I + 3]) + fabs(DX[I + 4]) + fabs(DX[I + 5]);
	}
	DASUM = DTEMP;
	return;

}

/*******************
 DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
 *     .. Scalar Arguments ..
 INTEGER INCX,N
 *     ..
 *     .. Array Arguments ..
 DOUBLE PRECISION DX(*)
 *     ..
 *
 *  Purpose
 *  =======
 *
 *     takes the sum of the absolute values.
 *     jack dongarra, linpack, 3/11/78.
 *     modified 3/93 to return if incx .le. 0.
 *     modified 12/3/93, array(1) declarations changed to array(*)
 *
 *
 *     .. Local Scalars ..
 DOUBLE PRECISION DTEMP
 INTEGER I,M,MP1,NINCX
 *     ..
 *     .. Intrinsic Functions ..
 INTRINSIC DABS,MOD
 *     ..
 DASUM = 0.0d0
 DTEMP = 0.0d0
 IF (N.LE.0 .OR. INCX.LE.0) RETURN
 IF (INCX.EQ.1) GO TO 20
 *
 *        code for increment not equal to 1
 *
 NINCX = N*INCX
 DO 10 I = 1,NINCX,INCX
 DTEMP = DTEMP + DABS(DX(I))
 10 CONTINUE
 DASUM = DTEMP
 RETURN
 *
 *        code for increment equal to 1
 *
 *
 *        clean-up loop
 *
 20 M = MOD(N,6)
 IF (M.EQ.0) GO TO 40
 DO 30 I = 1,M
 DTEMP = DTEMP + DABS(DX(I))
 30 CONTINUE
 IF (N.LT.6) GO TO 60
 40 MP1 = M + 1
 DO 50 I = MP1,N,6
 DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I+1)) + DABS(DX(I+2)) +
 +            DABS(DX(I+3)) + DABS(DX(I+4)) + DABS(DX(I+5))
 50 CONTINUE
 60 DASUM = DTEMP
 RETURN
 END

 ************/
