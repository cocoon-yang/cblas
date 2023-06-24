# include <stdio.h>
# include <stdlib.h>
# include <time.h>
#include <cfloat>
#include <cmath>
#include "cblas.h"


//
// https://netlib.org/lapack/explore-html/d8/dfd/dlas2_8f_source.html#l00106
//
//   \par Purpose:
//  =============
//   
//   DLAS2  computes the singular values of the 2-by-2 matrix
//      [  F   G  ]
//      [  0   H  ].
//   On return, SSMIN is the smaller singular value and SSMAX is the
//   larger singular value.
//  
//    Arguments:
//    ==========
//  
//   \param[in] F: double,  
//            The (0,0) element of the 2-by-2 matrix.
//
//   \param[in] G: double
//            The (0,1) element of the 2-by-2 matrix.
//  
//   \param[in] H: double, 
//            The (1,1) element of the 2-by-2 matrix.
//  
//   \param[out] SSMIN: double, 
//            The smaller singular value.
//  
//   \param[out] SSMAX: double, 
//            The larger singular value.
//  
//    Authors:
//    ========
//  
//   \author Univ. of Tennessee
//   \author Univ. of California Berkeley
//   \author Univ. of Colorado Denver
//   \author NAG Ltd.
//  
//   \ingroup OTHERauxiliary
//  
//   \par Further Details:
//    =====================
//  
//   \verbatim
//  
//    Barring over/underflow, all output quantities are correct to within
//    a few units in the last place (ulps), even in the absence of a guard
//    digit in addition/subtraction.
//  
//    In IEEE arithmetic, the code works correctly if ONE matrix element is
//    infinite.
//  
//    Overflow will not occur unless the largest singular value itself
//    overflows, or is within a few ulps of overflow. (On machines with
//    partial overflow, like the Cray, overflow may occur if the largest
//    singular value is within a factor of 2 of overflow.)
//  
//    Underflow is harmless if underflow is gradual. Otherwise, results
//    may correspond to a matrix modified by perturbations of size near
//    the underflow threshold.
//   \endverbatim
//  
//    =====================================================================
void dlas2(double f, double g, double h, double& ssmin, double& ssMAX)
{
	//  
	//    -- LAPACK auxiliary routine --
	//    -- LAPACK is a software package provided by Univ. of Tennessee,    --
	//    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
	//  
	//       .. Scalar Arguments ..
	//DOUBLE PRECISION   F, G, H, SSMAX, SSMIN
	//       ..
	//  
	//    ====================================================================
	//  
	//       .. Parameters ..
	double ZERO = 0.0;
	double ONE = 1.0;
	double TWO = 2.0;
	//       ..
	//       .. Local Scalars ..
	double as, at, au, c, fa, fhmn, fhmx, ga, ha;
	//  
	fa = fabs(f);
	ga = fabs(g);
	ha = fabs(h);
	fhmn = min(fa, ha);
	fhmx = MAX(fa, ha);

	if (fhmn == ZERO) {
		ssmin = ZERO;
		if (fhmx == ZERO) {
			ssMAX = ga;
		}
		else {
			ssMAX = MAX(fhmx, ga)*sqrt(ONE + (min(fhmx, ga) / (MAX(fhmx, ga))*MAX(fhmx, ga)));
		}
	}
	else {
		if (ga < fhmx) {
			as = ONE + fhmn / fhmx;
			at = (fhmx - fhmn) / fhmx;
			au = (ga / fhmx)*(ga / fhmx);
			c = TWO / (sqrt(as*as + au) + sqrt(at*at + au));
			ssmin = fhmn*c;
			ssMAX = fhmx / c;
		}
		else {
			au = fhmx / ga;

			if (au == ZERO) {
				//  
				// Avoid possible harmful underflow if expONEnt range
				// asymmetric (true SSMIN may not underflow even if
				// AU underflows)
				//  
				ssmin = (fhmn*fhmx) / ga;
				ssMAX = ga;
			}
			else {
				as = ONE + fhmn / fhmx;
				at = (fhmx - fhmn) / fhmx;
				c = ONE / (sqrt(ONE + (as*au)*(as*au)) + sqrt(ONE + (at*au)*(at*au)));
				ssmin = (fhmn*c)*au;
				ssmin = ssmin + ssmin;
				ssMAX = ga / (c + c);
			}
		}
	}
	return;
}

//
// DPSTF2 computes the Cholesky factorization with complete
// pivoting of a real symmetric positive semidefinite matrix A.
//
// The factorization has the form
//    P**T * A * P = U**T * U ,  if UPLO = 'U',
//    P**T * A * P = L  * L**T,  if UPLO = 'L',
// where U is an upper triangular matrix and L is lower triangular, and
// P is stored as vector PIV.
//
// This algorithm does not attempt to check that A is positive
// semidefinite. This version of the algorithm calls level 2 BLAS.
//
//          UPLO is CHARACTER*1
//          Specifies whether the upper or lower triangular part of the
//          symmetric matrix A is stored.
//          = 'U':  Upper triangular
//          = 'L':  Lower triangular

//          N is INTEGER
//          The order of the matrix A.  N >= 0.

//          A is DOUBLE PRECISION array, dimension (LDA,N)
//          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
//          n by n upper triangular part of A contains the upper
//          triangular part of the matrix A, and the strictly lower
//          triangular part of A is not referenced.  If UPLO = 'L', the
//          leading n by n lower triangular part of A contains the lower
//          triangular part of the matrix A, and the strictly upper
//          triangular part of A is not referenced.
//
//          On exit, if INFO = 0, the factor U or L from the Cholesky
//          factorization as above.

//          PIV is INTEGER array, dimension (N)
//          PIV is such that the nonzero entries are P( PIV(K), K ) = 1.

//          RANK is INTEGER
//          The rank of A given by the number of steps the algorithm
//          completed.

//          TOL is DOUBLE PRECISION
//          User defined tolerance. If TOL < 0, then N*U*MAX( A( K,K ) )
//          will be used. The algorithm terminates at the (K-1)st step
//          if the pivot <= TOL.

//          WORK is DOUBLE PRECISION array, dimension (2*N)
//          Work space.

//          INFO is INTEGER
//          < 0: If INFO = -K, the K-th argument had an illegal value,
//          = 0: algorithm completed successfully, and
//          > 0: the matrix A is either rank deficient with computed rank
//               as returned in RANK, or is not positive semidefinite. See
//               Section 7 of LAPACK Working Note #161 for further
//               information.
void dpstf2(char UPLO, int N, double* A, int LDA, int* PIV, int RANK, double TOL, double* WORK, int INFO)
{
	// Local parameters
	double ZERO = 0.0;

	INFO = 0;
	bool UPPER = cblas_lsame(UPLO, 'U');
	if (!UPPER && !cblas_lsame(UPLO, 'L')) {
		INFO = -1;
	}
	else if (N < 0) {
		INFO = -2;
	}
	else if (LDA < MAX(1, N)) {
		INFO = -4;
	}
	if (INFO != 0) {
		xerbla("DPSTF2", -INFO);
		return;
	}

	// Quick return 
	if (N == 0)
	{
		return;
	}

	// Initialize PIV
	for (int i = 0; i < N; i++)
	{
		PIV[i] = i;
	}
 
	// Compute stopping value
	//DO I = 2, N
	//IF(A(I, I).GT.AJJ) THEN
	//PVT = I
	//AJJ = A(PVT, PVT)
	//END IF
	//END DO
	int PVT = 1;
	double AJJ = A[PVT * LDA + PVT];
	for (int i = 1; i < N; i++)
	{
		if (A[PVT * LDA + PVT] > AJJ)
		{
			PVT = i;
			AJJ = A[PVT* LDA + PVT];
		}
	}

	if ((AJJ < ZERO) || std::isnan(AJJ)) {
		RANK = 0;
		INFO = 1;
		return;
	}
	// to be continue
}