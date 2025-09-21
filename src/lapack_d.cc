# include <stdio.h>
# include <stdlib.h>
# include <time.h>
#include <cfloat>
#include <cmath>
#include "cblas.h"


//#define min(x,y) (((x) < (y)) ? (x) : (y))
//#define max(x,y) (((x) > (y)) ? (x) : (y)) 

extern double sign(double A, double B);
extern double dlamch(const char* CMACH); 
extern double dlapy2(double x, double y); 
extern int iladlc(int m, int n, double* pA, int lda);
extern int iladlr(int m, int n, double* pA, int lda); 
extern int ilaenv(int ispec, const char* name, char OPTS, int N1, int N2, int N3, int N4);
extern double cblas_idamax(int n, double* dx, int incx);
extern void cblas_xerbla(const char* SRNAME, int INFO);
extern bool cblas_lsame(char CA, char CB);
extern void cblas_daxpy(const int N, const double alpha, const double* X,
	const int incX, double* Y, const int incY);
extern void cblas_dscal(const int N, const double DA, double* DX, const int INCX);
extern int cblas_dgemv(char TRANS, int M, int N, double ALPHA, double* A, int LDA,
	double* X, int INCX, double BETA, double* Y, int INCY);
extern int cblas_dger(int M, int N, double ALPHA, double* X, int INCX, double* Y,
	int INCY, double* A, int LDA);
extern int iladlr(int m, int n, double* pA, int lda); 
extern int cblas_iladlr(int m, int n, double* pA, int lda);
extern double cblas_dnrm2(const int N, const double* X, const int incX);
extern void cblas_dtrsm(char side, char uplo, char transa, char diag, int m, int n,
	double alpha, double a[], int lda, double b[], int ldb);
extern void cblas_dlaswp(int n, double* a, int lda, int k1, int k2, int* ipiv, int incx);
extern void cblas_dgemm(char transa, char transb, int m, int n, int k,
	double alpha, double a[], int lda, double b[], int ldb, double beta,
	double c[], int ldc);


/***
=========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE DGETRF2( M, N, A, LDA, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       int            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       int            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * )
*       ..
*
*
@brief Purpose:
=============

\verbatim

DGETRF2 computes an LU factorization of a general M-by-N matrix A
using partial pivoting with row interchanges.

The factorization has the form
A = P * L * U
where P is a permutation matrix, L is lower triangular with unit
diagonal elements (lower trapezoidal if m > n), and U is upper
triangular (upper trapezoidal if m < n).

This is the recursive version of the algorithm. It divides
the matrix into four submatrices:

[  A11 | A12  ]
[ -----|----- ] = A
[  A21 | A22  ]
where A11 is n1 by n1 and A22 is n2 by n2
with n1 = min(m,n)/2, n2 = n-n1

The subroutine calls itself to factor
[ A11 ]
[ --- ],
[ A12 ]

do the swaps on
[ A12 ]
[ --- ], solve A12, update A22,
[ A22 ]

then calls itself to factor A22 and do the swaps on A21.

\endverbatim
*
*  Arguments:
*  ==========
*
\param[in] M: int
The number of rows of the matrix A.  M >= 0.

\param[in] N: int
The number of columns of the matrix A.  N >= 0.

\param[in,out] A
A is double array, dimension (LDA,N)
On entry, the M-by-N matrix to be factored.
On exit, the factors L and U from the factorization
A = P*L*U; the unit diagonal elements of L are not stored.

\param[in] LDA: int
The leading dimension of the array A.  LDA >= max(1,M).

\param[out] IPIV: IPIV is int array, dimension (min(M,N))
The pivot indices; for 1 <= i <= min(M,N), row i of the
matrix was interchanged with row IPIV(i).

\param[out] INFO: int
= 0:  successful exit
< 0:  if INFO = -i, the i-th argument had an illegal value
> 0:  if INFO = i, U(i,i) is exactly zero. The factorization
has been completed, but the factor U is exactly
singular, and division by zero will occur if it is used
to solve a system of equations.
*
*  Authors:
*  ========
*
\author Univ. of Tennessee
\author Univ. of California Berkeley
\author Univ. of Colorado Denver
\author NAG Ltd.
*
\ingroup doubleGEcomputational
*
*  =====================================================================
RECURSIVE SUBROUTINE dgetrf2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  =====================================================================
********/
void dgetrf2(int m, int n, double* pA, int LDA, int* ipiv, int& INFO)
{
	// Local variables
	double ZERO = 0.0;
	double ONE = 1.0;
	int  n1, n2, i;
	int iinfo = 0;
	double sfmin, temp;

	// Test the input parameters
	INFO = 0;
	if (m < 0) {
		INFO = -1;
	}
	else if (n < 0) {
		INFO = -2;
	}
	else if (nullptr == pA) {
		INFO = -3;
	}
	else if (LDA < std::max(1, n)) {
		INFO = -4;
	}
	else if (nullptr == ipiv) {
		INFO = -5;
	}
	if (INFO != 0) {
		// Show the error code
		cblas_xerbla("DGETRF2", -INFO);
		return;
	}

	//Quick return if possible 
	if (m == 0 || n == 0)
	{
		return;
	}

	MData A(m, n, pA, LDA);
	//A.setData(pA);

	if (m == 1) {
		//  Use unblocked code for one row case
		//  Just need to handle IPIV and INFO 
		ipiv[0] = 0;
		if (A[0][0] == ZERO) {
			INFO = 1;
		}
	}
	else if (n == 1) {
		//  Use unblocked code for one column case 

		// Compute machine safe minimum 
		sfmin = dlamch("S");
		// Find pivot and test for singularity 
		// for the i-th column 
		i = cblas_idamax(m, pA, LDA);
		ipiv[0] = i;
		if (A[i][0] != ZERO) {
			// Apply the interchange 
			if (i != 0) {
				temp = A[0][0];
				A[0][0] = A[i][0];
				A[i][0] = temp;
			}

			// Compute elements 1:M-1 of the column 
			if (fabs(A[0][0]) >= sfmin) {
				//cblas_dscal(m - 1, ONE / A[0][0], pA + (1 * LDA), LDA); 
				cblas_dscal(m - C_VEC_OFFSET, ONE / A[0][0], A.sub(1, 0), LDA);
			}
			else {
				for (int i = 1; i < m; i++)
				{
					A[i][0] = A[i][0] / A[0][0];
				}
			}
		}
		else {
			INFO = 1;
		}
	}
	else {
		// Block code 
		//  Use recursive code 
		n1 = std::min(m, n) / 2;
		n2 = n - n1;
		//
		//        [ A11 ]
		// Factor [ --- ]
		//        [ A21 ]
		//
		dgetrf2(m, n1, pA, LDA, ipiv, iinfo);

		if (INFO == 0 & iinfo > 0)
		{
			// Something goes wrong in the Block factorization 
			INFO = iinfo;
		}
		//
		//                              [ A12 ]
		//        Apply interchanges to [ --- ]
		//                              [ A22 ]
		//
		cblas_dlaswp(n2, A.sub(0, n1), LDA, 0, n1 - C_VEC_OFFSET, ipiv, 1);
		//
		//        Solve A12
		//
		cblas_dtrsm('L', 'L', 'N', 'U', n1, n2, ONE, pA, LDA, A.sub(0, n1), LDA);
		//
		//        Update A22
		//
		cblas_dgemm('N', 'N', m - n1, n2, n1, -ONE, A.sub(n1, 0), LDA, A.sub(0, n1), LDA, ONE, A.sub(n1, n1), LDA);
		//
		//        Factor A22
		//
		dgetrf2(m - n1, n2, A.sub(n1, n1), LDA, ipiv + n1, iinfo);
		//
		//        Adjust INFO and the pivot indices
		//
		if ((INFO == 0) && (iinfo > 0))
		{
			INFO = iinfo + n1;
		}
		for (int i = n1; i < std::min(m, n); i++) {
			ipiv[i] = ipiv[i] + n1;
		}
		//
		//        Apply interchanges to A21
		//
		cblas_dlaswp(n1, pA, LDA, n1, std::min(m, n) - C_VEC_OFFSET, ipiv, 1);
	}
	return;
}


void clapack_dgetrf(int m, int n, double* pA, int lda, int* ipiv, int& info)
{
	// Temperate variables 
	double one = 1.0;
	int i, iinfo, j, jb, nb;

	iinfo = 0;

	/*
	*Test the input parameters.
	*/
	info = 0;
	if (m < 0) {
		info = -1;
	}
	else if (n < 0) {
		info = -2;
	}
	else if (lda < std::max(1, m)) {
		info = -4;
	}
	if (info != 0) {
		//xerbla("DGETRF", -info);
		return;
	}

	/*
	*  Quick return if possible
	*/
	if ((m == 0) || (n == 0))
	{
		return;
	}

	MData a(m, lda, pA, lda);
	//a.setData(pA);


	/*
	*     Determine the block size for this environment.
	*/
	nb = ilaenv(1, "DGETRF", ' ', m, n, -1, -1);

	// 
	//printf("     \n");
	//std::cout << "dgetrf() block size: " << nb << std::endl;
	//



	if ((nb <= 1) || (nb >= std::min(m, n))) {
		/*
		* Use unblocked code.
		*/
		dgetrf2(m, n, pA, lda, ipiv, info);
	}
	else {
		/*
		*  Use blocked code.
		*/
		for (j = 0; j < std::min(m, n); j += nb)
		{
			jb = std::min(std::min(m, n) - j + 1, nb);

			/*
			*    Factor diagonal and subdiagonal blocks and test for exact
			*    singularity.
			*/
			dgetrf2(m - j, jb, a.sub(j, j), lda, ipiv + j, iinfo);

			/*
			*   Adjust INFO and the pivot indices.
			*/
			if ((info == 0) && (iinfo > 0))
			{
				info = iinfo + j - 1;
			}
			for (i = j; i < std::min(m, j + jb - 1); i++)
			{
				ipiv[i] = j - 1 + ipiv[i];
			}

			/*
			*   Apply interchanges to columns 1:J-1.
			*/
			cblas_dlaswp(j - 1, pA, lda, j, j + jb - 1, ipiv, 1);

			if ((j + jb) <= n)
			{
				/*
				*   Apply interchanges to columns J+JB:N.
				*/
				cblas_dlaswp(n - j - jb + 1, a.sub(0, j + jb), lda, j, j + jb - 1, ipiv, 1);

				/*
				*   Compute block row of U.
				*/
				cblas_dtrsm('L', 'L', 'N', 'U', jb,
					n - j - jb + 1, one, a.sub(j, j), lda, a.sub(j, j + jb),
					lda);

				if ((j + jb) <= m)
				{
					/*
					*  Update trailing submatrix.
					*/
					cblas_dgemm('N', 'N', m - j - jb + 1,
						n - j - jb + 1, jb, -one, a.sub(j + jb, j), lda,
						a.sub(j, j + jb), lda, one, a.sub(j + jb, j + jb),
						lda);
				}
			}
		}
	}
}



void clapack_dgetrs(char trans, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb, int info)
{
	double one = 1.0;
	// ..Local Scalars ..
	bool notran;

	// Test the input parameters.
	info = 0;
	notran = cblas_lsame(trans, 'N');
	if ((!notran) & (!cblas_lsame(trans, 'T')) & (!cblas_lsame(trans, 'C'))) {
		info = -1;
	}
	else if (n < 0) {
		info = -2;
	}
	else if (nrhs < 0) {
		info = -3;
	}
	//else if (lda < std::max(1, n)) {
	//	info = -5;
	//}
	//else if (ldb < std::max(1, n)) {
	//	info = -8;
	//}

	//printf(" cblas_dgetrs:   n = %i  \n", n);
	//printf(" cblas_dgetrs: ldb = %i  \n", ldb);


	if (info != 0) {
		cblas_xerbla("DGETRS", -info);
		return;
	}

	/*
	*  Quick return if possible
	*/
	if ((n == 0) || (nrhs == 0)) {
		return;
	}

	if (notran) {


		//printf(" cblas_dgetrs: notran  \n");

		/*
		* Solve A* X = B.
		*
		*Apply row interchanges to the right hand sides.
		*/
		cblas_dlaswp(nrhs, b, ldb, 0, n - 1, ipiv, 1);

		/*
		*Solve L * X = B, overwriting B with X.
		*/
		cblas_dtrsm('L', 'L', 'N', 'U', n,
			nrhs, one, a, lda, b, ldb);



		/*
		*Solve U * X = B, overwriting B with X.
		*/
		cblas_dtrsm('L', 'U', 'N', 'N', n,
			nrhs, one, a, lda, b, ldb);


		//showMatrix_d(b, n, 1);

	}
	else {
		/*
		*Solve A * *T * X = B.
		*
		*Solve U * *T * X = B, overwriting B with X.
		*/
		cblas_dtrsm('L', 'U', 'T', 'N', n,
			nrhs, one, a, lda, b, ldb);
		/*
		*Solve L * *T * X = B, overwriting B with X.
		*/
		cblas_dtrsm('L', 'L', 'T', 'U', n, nrhs, one, a, lda, b, ldb);
		/*
		*Apply row interchanges to the solution vectors.
		*/
		cblas_dlaswp(nrhs, b, ldb, 0, n - 1, ipiv, -1);
	}
	return;
}

///
// DGESV computes the solution to a real system of linear equations
// A* X = B,
// where A is an N - by - N matrix and X and B are N - by - NRHS matrices.
//
// The LU decomposition with partial pivoting and row interchanges is
// used to factor A as
// A = P * L * U,
// where P is a permutation matrix, L is unit lower triangular, and U is
// upper triangular.The factored form of A is then used to solve the
// system of equations A* X = B.
//
void clapack_dgesv(int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb, int info)
{
	int INFO, LDA, LDB, N, NRHS;
	//
	//Array Arguments ..
	int* IPIV;
	double* A;
	double* B;
	/* ..
	*
	* ====================================================================
	*
	* Test the input parameters.
	*/
	info = 0;
	if (n < 0) {
		info = -1;
	}
	else if (nrhs < 0) {
		info = -2;
	}
	//else if (lda < std::max(1, n)) {
	//	info = -4;
	//}
	//else if (ldb < std::max(1, n)) {
	//	info = -7;
	//}


	printf("DGESV: info: %i \n", info);


	if (info != 0) {
		cblas_xerbla("DGESV", -info);
		return;
	}
	/*
	*Compute the LU factorization of A.
	*/
	clapack_dgetrf(n, n, a, lda, ipiv, info);
	if (info == 0) {
		/*
		*Solve the system A * X = B, overwriting B with X.
		*/
		clapack_dgetrs('N', n, nrhs, a, lda, ipiv, b, ldb, info);
	}
	return;
}



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
	fhmn = std::min(fa, ha);
	fhmx = std::max(fa, ha);

	if (fhmn == ZERO) {
		ssmin = ZERO;
		if (fhmx == ZERO) {
			ssMAX = ga;
		}
		else {
			ssMAX = std::max(fhmx, ga)*sqrt(ONE + (std::min(fhmx, ga) / (std::max(fhmx, ga))* std::max(fhmx, ga)));
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
	else if (LDA < std::max(1, N)) {
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



/****
\brief DLARF1F applies a real elementary reflector H to a real m by n matrix
 C, from either the left or the right. H is represented in the form

 H = I - tau * v * v * *T

 where tau is a real scalar and v is a real vector.

 If tau = 0, then H is taken to be the unit matrix.

@param[in] side, char
		= 'L': form  H * C
		= 'R': form  C * H
@param[in] m, int  The number of rows of the matrix C.
@param[in] n, int, The number of columns of the matrix C.
@param[in] v, double*,  dimension
					 (1 + (M-1)*abs(INCV)) if SIDE = 'L'
				  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
		  The vector v in the representation of H. V is not used if
		  TAU = 0. V(1) is not referenced or modified.
@param[in] incv, int, The increment between elements of v. INCV <> 0.
@param[in] tau, double, The value tau in the representation of H.
@param[in] c, double*, dimension (LDC,N)
		  On entry, the m by n matrix C.
		  On exit, C is overwritten by the matrix H * C if SIDE = 'L',
		  or C * H if SIDE = 'R'.
@param[in] ldc,  int, The leading dimension of the array C. LDC >= max(1,M).
@param[in] work, double*, dimension
						 (N) if SIDE = 'L'
					  or (M) if SIDE = 'R'
*/
void clapack_dlarf1f(char side, int m, int n, double* v,
	int incv, double tau, double* c, int ldc, double* work)
{
	double one = 1.0;
	double zero = 0.0;

	bool applyleft = cblas_lsame(side, 'L');
	int lastv = 1;
	int lastc = 0;
	int i = 0;
	if (tau != zero) {
		// Set up variables for scanning V.
		// LASTV begins pointing to the end of V.
		if (applyleft) {
			lastv = m;
		}
		else {
			lastv = n;
		}
		if (incv > 0) {
			i = 1 + (lastv - 1) * incv;
		}
		else {
			i = 1;
		}
		//Look for the last non - zero row in V.
		//Since we are assuming that V[0] = 1, and it is not stored, so we
		//shouldn't access it.
		while ((lastv > 1) & (v[i] == zero)) {
			lastv = lastv - 1;
			i = i - incv;
		}
		if (applyleft) {
			//Scan for the last non - zero column in C(1:lastv, : ).
			lastc = cblas_iladlr(lastv, n, c, ldc);
			lastc = lastc - 1;
		}
		else {
			// Scan for the last non - zero row in C(:, 1 : lastv).
			lastc = iladlr(m, lastv, c, ldc);
			//lastc = lastc;
		}
	}
	if (lastc == 0) {
		return;
	}
	if (applyleft) {
		/*
		* Form  H * C
		*/
		//Check if lastv = 1. This means v = 1, So we just need to compute
		// C : = HC = (1 - \tau)C.
		if (lastv == 1) {
			/*
			* C(1, 1:lastc) : = (1 - tau) * C(1, 1:lastc)
			*/
			cblas_dscal(lastc, one - tau, c, ldc);
		}
		else {
			/*
			* w(1:lastc, 1) : = C(1:lastv, 1 : lastc) * *T * v(1:lastv, 1)
			*/
			// w(1:lastc, 1) : = C(2:lastv, 1 : lastc) * *T * v(2:lastv, 1)
			cblas_dgemv('T', lastv - 1, lastc, one, c + (ldc),
				ldc, v + incv, incv, zero, work, 1);

			//
			//std::cout << " C:" << std::endl;
			//showMatrix(c, m, n);
			//std::cout << "work:" << std::endl;
			//showVector(work, m);
			//


			//w(1:lastc, 1) += C(1, 1:lastc) * *T * v(1, 1) = C(1, 1:lastc) * *T
			cblas_daxpy(lastc, one, c, 1, work, 1);
			// 
			//std::cout << " C:" << std::endl;
			//showMatrix(c, m, n);
			//std::cout << "work:" << std::endl;
			//showVector(work, m);
			//

			/*
			* C(1:lastv, 1 : lastc) := C(...) - tau * v(1:lastv, 1) * w(1:lastc, 1) * *T
			*/
			// C(1, 1:lastc) : = C(...) - tau * v(1, 1) * w(1:lastc, 1) * *T
			// = C(...) - tau * w(1:lastc, 1) * *T
			cblas_daxpy(lastc, -tau, work, 1, c, 1);
			// 
			//std::cout << " C:" << std::endl;
			//showMatrix(c, m, n);
			//std::cout << "work:" << std::endl;
			//showVector(work, m);
			//std::cout << "v:" << std::endl;
			//showVector(v + incv, m);
			//

			//
			//double tmp[]{
			//	0, 0, 1, 1
			//};
			//cblas_dger(lastv - 1, lastc, -tau, v + incv, incv, work, 1,
			//	tmp, 2);


			//std::cout << " tmp:" << std::endl;
			//showMatrix(tmp, 2, 2);
			//


			//C(2:lastv, 1 : lastc) : = C(...) - tau * v(2:lastv, 1) * w(1:lastc, 1) * *T
			cblas_dger(lastv - 1, lastc, -tau, v + incv, incv, work, 1,
				c + (ldc), ldc);
			// 
			//std::cout << " C:" << std::endl;
			//showMatrix(c, m, n);
			// 
		}
	}
	else {
		/*
		*Form  C * H
		*/
		// Check if n = 1. This means v = 1, so we just need to compute
		// C : = CH = C(1 - \tau).
		if (lastv == 1) {
			/*
			*C(1:lastc, 1) : = (1 - tau) * C(1:lastc, 1)
			*/
			cblas_dscal(lastc, one - tau, c, 1);
		}
		else {
			/*
			*w(1:lastc, 1) : = C(1:lastc, 1 : lastv) * v(1:lastv, 1)
			*/
			//w(1:lastc, 1) : = C(1:lastc, 2 : lastv) * v(2:lastv, 1)
			cblas_dgemv('N', lastv - 1, lastc, one,  // 
				c + 1, ldc, v + (1), incv, zero, work, 1);

			//
			//std::cout << "work:" << std::endl;
			//showVector(work, m);
			//


			//w(1:lastc, 1) += C(1:lastc, 1) v(1, 1) = C(1:lastc, 1)
			cblas_daxpy(lastc, one, c, ldc, work, 1);

			//
			//std::cout << "work:" << std::endl;
			//showVector(work, m); 
			//std::cout << " C:" << std::endl;
			//showMatrix(c, m, n);
			// 



			/*
			*C(1:lastc, 1 : lastv) := C(...) - tau * w(1:lastc, 1) * v(1:lastv, 1) * *T
			*/
			//C(1:lastc, 1) : = C(...) - tau * w(1:lastc, 1) * v(1, 1) * *T
			//= C(...) - tau * w(1:lastc, 1)
			cblas_daxpy(lastc, -tau, work, 1, c, ldc);

			// 
			//std::cout << " C:" << std::endl;
			//showMatrix(c, m, n);
			// 

			//C(1:lastc, 2 : lastv) : = C(...) - tau * w(1:lastc, 1) * v(2:lastv) * *T
			cblas_dger(lastc, lastv, -tau, work, 1, v + (incv),
				incv, c + (1), ldc);
		}
	}
	return;
}


// 
// QR Fraction 
//


/***
* \brief DLARFG generates a real elementary reflector H of order n, such
* that
*
* H* (alpha) = (beta), H** T* H = I.
*      (x)       (0)
*
* where alpha and beta are scalars, and x is an(n - 1) - element real
* vector.H is represented in the form
*
* H = I - tau * (1) * (1 v * *T),
*               (v)
*
* where tau is a real scalar and v is a real(n - 1) - element
* vector.
*
* If the elements of x are all zero, then tau = 0 and H is taken to be
* the unit matrix.
*
* Otherwise  1 <= tau <= 2.
*
*  @param[in] n: int, The order of the elementary reflector.
*  @param[in/out] alpha: double&, On entry, the value alpha.
*                            On exit, it is overwritten with the value beta.
*  @param[in] x: double*, On entry, the vector x.
*                         On exit, it is overwritten with the vector v.
*  @param[in] incx: int, The increment between elements of X. INCX > 0.
*  @param[out] tau: double, The value tau.
*
*/
void clapack_dlarfg(int	n, double& alpha, double* x, int incx, double& tau)
{
	double zero = 0.0;
	double one = 1.0;

	if (n <= 1) {
		tau = zero;
		return;
	}

	double xnorm = cblas_dnrm2(n - 1, x, incx);

	if (xnorm == zero) {
		/*
		* H = I
		*/
		tau = zero;
	}
	else {
		/*
		*general case
		*/
		// beta = sqrt(alpha^2 + xnorm^2 )= norm( x )
		double beta = -sign(dlapy2(alpha, xnorm), alpha);
		double safmin = dlamch("S") / dlamch("E");
		int knt = 0;
		if (abs(beta) < safmin) {
			/*
			* XNORM, BETA may be inaccurate; scale X and recompute them
			*/
			double rsafmn = one / safmin;
			//10       CONTINUE 
			do {
				knt = knt + 1;
				cblas_dscal(n - 1, rsafmn, x, incx);
				beta = beta * rsafmn;
				alpha = alpha * rsafmn;
			} while ((fabs(beta) < safmin) & (knt < 20));
			//$         GO TO 10
			/*
			*  New BETA is at most 1, at least SAFMIN
			*/
			xnorm = cblas_dnrm2(n - 1, x, incx);
			beta = -sign(dlapy2(alpha, xnorm), alpha);
		}
		tau = (beta - alpha) / beta;
		cblas_dscal(n - 1, one / (alpha - beta), x, incx);
		/*
		*If ALPHA is subnormal, it may lose relative accuracy
		*/
		for (int j = 0; j < knt; j++) {
			beta = beta * safmin;
			// 20      CONTINUE
		}
		alpha = beta;
	}

	return;
}

/**
 * @brief  DLARF applies a real elementary reflector H to a real m by n matrix
 *  C, from either the left or the right. H is represented in the form
 *
 *       H = I - tau * v * v**T
 *
 * where tau is a real scalar and v is a real vector.
 *
 * If tau = 0, then H is taken to be the unit matrix.
 *
 *
 * @param[in] side: char,  = 'L': form  H * C
 *                     = 'R': form  C * H
 * @param[in] m: int, The number of rows of the matrix C.
 * @param[in] n: int, The number of columns of the matrix C.
 * @param[in] v: double array, dimension
 *                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
 *                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
 *          The vector v in the representation of H. V is not used if
 *          TAU = 0.
 * @param[in] incv: int, The increment between elements of v. INCV != 0.
 * @param[in] tau: double, The value tau in the representation of H.
 * @param C: double array, dimension (LDC,N)
 *          On entry, the m by n matrix C.
 *          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
 *          or C * H if SIDE = 'R'.
 * @param[in] ldc: int, The leading dimension of the array C. LDC >= max(1,M).
 * @param[out] work: double array, dimension
 *                         (N) if SIDE = 'L'
 *                      or (M) if SIDE = 'R'
 *
*/
void clapack_dlarf(char side, int m, int n, double* v, int incv,
	double tau, double* c, int ldc, double* work)
{
	// Temperate variables
	double one = 1.0;
	double zero = 0.0;
	int i = 0;

	bool applyleft = cblas_lsame(side, 'L');
	int lastv = 0;
	int	lastc = 0;

	if (tau != zero) {
		//Set up variables for scanning V.LASTV begins pointing to the end
		//of V.
		if (applyleft) {
			lastv = m;
		}
		else {
			lastv = n;
		}
		if (incv > 0) {
			i = 1 + (lastv - 1) * incv;
		}
		else {
			i = 1;
		}
		//Look for the last non-zero row in V.
		while ((lastv > 0) & (v[i] == zero)) {
			lastv = lastv - 1;
			i = i - incv;
		}
		if (applyleft) {
			//Scan for the last non - zero column in C[0:lastv, : ].
			lastc = iladlc(lastv, n, c, ldc);
		}
		else {
			//Scan for the last non - zero row in C[:, 0 : lastv].
			lastc = iladlr(m, lastv, c, ldc);
		}
	}

	// Note that lastc.eq.0 renders the BLAS operations null; no special
	// case is needed at this level.
	if (applyleft) {
		/*
		* Form  H * C
		*/
		if (lastv > 0) {
			/*
			* w(1:lastc, 1) := C(1:lastv, 1 : lastc) * *T * v(1:lastv, 1)
			*/
			cblas_dgemv('T', lastv, lastc, one, c, ldc, v,
				incv, zero, work, 1);
			/*
			*C(1:lastv, 1 : lastc) := C(...) - v(1:lastv, 1) * w(1:lastc, 1) * *T
			*/
			cblas_dger(lastv, lastc, -tau, v, incv, work, 1, c, ldc);
		}
	}
	else {
		/*
		*Form  C * H
		*/
		if (lastv > 0) {
			/*
			*w(1:lastc, 1) := C(1:lastc, 1 : lastv) * v(1:lastv, 1)
			*/
			cblas_dgemv('N', lastc, lastv, one, c, ldc,
				v, incv, zero, work, 1);
			/*
			*C(1:lastc, 1 : lastv) : = C(...) - w(1:lastc, 1) * v(1:lastv, 1) * *T
			*/
			cblas_dger(lastc, lastv, -tau, work, 1, v, incv, c, ldc);
		}
	}
} 

/****
 \brief DGEQR2 computes a QR factorization of a real m - by - n matrix A :

 A = Q * (R),
		 (0)

 where:

 Q is a m - by - m orthogonal matrix;
 R is an upper - triangular n - by - n matrix;
 0 is a(m - n) - by - n zero matrix, if m > n.

  The matrix Q is represented as a product of elementary reflectors

	 Q = H(1) H(2) . . . H(k), where k = min(m,n).

  Each H(i) has the form

	 H(i) = I - tau * v * v**T

  where tau is a real scalar, and v is a real vector with
  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
  and tau in TAU(i).

  @param[in] m: int, The number of rows of the matrix A.  M >= 0.

  @param[in] n: int, The number of columns of the matrix A.  N >= 0.

  @param[out] A: double array, dimension (LDA,N)
		  On entry, the m by n matrix A.
		  On exit, the elements on and above the diagonal of the array
		  contain the min(m,n) by n upper trapezoidal matrix R (R is
		  upper triangular if m >= n); the elements below the diagonal,
		  with the array TAU, represent the orthogonal matrix Q as a
		  product of elementary reflectors (see Further Details).

  @param[in] lda: int, The leading dimension of the array A.  LDA >= max(1,M).

  @param[out] work: double array, dimension (min(M,N))
		   The scalar factors of the elementary reflectors (see Further
		   Details).

  @param[out] work: double array, dimension (N)

  @param[out] info
		   = 0: successful exit
		   < 0: if INFO = -i, the i-th argument had an illegal value
 */

void clapack_dgeqr2(int m, int n, double* a, int lda, double* tau, double* work, int& info)
{
	double one = 1.0;
	info = 0;
	if (m < 0) {
		info = -1;
	}
	else if (n < 0) {
		info = -2;
	}
	else if (lda > (std::max)(1, m)) {
		info = -4;
	}
	if (info != 0) {
		cblas_xerbla("DGEQR2", -info);
		return;
	}

	MData A(m, n, a, lda);
	//A.setData(a);
	int k = (std::min)(m, n);

	for (int i = 0; i < k - 1; i++) {
		/*
		*Generate elementary reflector H(i) to annihilate A(i + 1:m, i)
		*/
		int sub_m = (std::min)(i, m - 1);
		clapack_dlarfg(m - i, A[i][i], A.sub(i + 1, sub_m), lda, tau[i]);


		//dlarfg(m, A[0], A + incA, incA, tau);

		//
		A.show();
		//

		if (i < (n - 1)) {
			/*
			*Apply H(i) to A(i:m, i + 1 : n) from the left
			*/
			clapack_dlarf1f('L', m - i, n - i - 1, A.sub(i, i), lda, tau[i], A.sub(i, i + 1), lda, work);


			// dlarf1f('L', m - i + 1, n - i, A, n, tau, A + 1, incA, work);

			//
			A.show();
			//
		}
	}
	return;
}
