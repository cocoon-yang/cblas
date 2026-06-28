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
extern  void cblas_dcopy(const int N, const double* X, const int incX, double* Y, const int incY);
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
extern void cblas_dtrmm(char side, char uplo, char transa, char diag, int m, int n,
	double alpha, double a[], int lda, double b[], int ldb);
extern void cblas_dsyrk(char uplo, char trans, int n, int k,
	double alpha, double* A, int lda, double beta, double* c, int ldc);
extern void cblas_dtrmv(char uplo, char trans, char diag, int n, double a[], int lda,
	double x[], int incx);

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



void clapack_dgetrs(char trans, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb, int& info)
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
void clapack_dgesv(int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb, int& info)
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

	//
	//printf("   \n");
	//printf("   dlarf1f()  \n");
	//printf("       tau: %5.3f \n", tau);
	//printf("    column: %i \n", m);
	//printf("       row: %i \n", n);
	//MData C(m, n, c, ldc);
	//C.show();
	//std::cout << "work:" << std::endl;
	//showVector(work, m);
	//

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
			//Scan for the last non - zero column in C(0:lastv, : ).
			//lastc = iladlr(lastv, n, c, 1); 

			lastc = iladlc(lastv, n, c, ldc);
			lastc = lastc;
		}
		else {
			// Scan for the last non - zero row in C(:, 0 : lastv).
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
			//C.show();
			//std::cout << "work:" << std::endl;
			//showVector(work, m);
			//


			//w(1:lastc, 1) += C(1, 1:lastc) * *T * v(1, 1) = C(1, 1:lastc) * *T
			cblas_daxpy(lastc, one, c, 1, work, 1);
			// 
			//std::cout << " C:" << std::endl;
			//C.show();
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
			//C.show();
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
			//C.show();
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
			cblas_dgemv('N', lastv, lastc, one,  // 
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
			//


			// 
			//std::cout << " C:" << std::endl;
			//C.show();
			// 



			/*
			*C(1:lastc, 1 : lastv) := C(...) - tau * w(1:lastc, 1) * v(1:lastv, 1) * *T
			*/
			//C(1:lastc, 1) : = C(...) - tau * w(1:lastc, 1) * v(1, 1) * *T
			//= C(...) - tau * w(1:lastc, 1)
			cblas_daxpy(lastc, -tau, work, 1, c, ldc);

			// 
			//std::cout << " C:" << std::endl;
			//C.show();
			// 

			//C(1:lastc, 2 : lastv) : = C(...) - tau * w(1:lastc, 1) * v(2:lastv) * *T
			cblas_dger(lastc, lastv, -tau, work, 1, v + (incv),
				incv, c + (1), ldc);
		}
	}
	return;
}


// 
// Cholesky Fraction 
//


/***
DPOTRF2 computes the Cholesky factorization of a real symmetric
positive definite matrix A using the recursive algorithm.

The factorization has the form
A = U^T U,  if UPLO = 'U', or
A = L L^T,  if UPLO = 'L',
where U is an upper triangular matrix and L is lower triangular.

This is the recursive version of the algorithm. It divides
the matrix into four submatrices:

	[  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
A = [ -----|----- ]  with n1 = n/2
	[  A21 | A22  ]       n2 = n-n1

The subroutine calls itself to factor A11. Update and scale A21
or A12, update A22 then calls itself to factor A22.
*/
void clapack_dpotrf2(char uplo, int n, double* A, int lda, int info)
{
	double one = 1.0;
	double zero = 0.0;

	int i, j, nb;
	int jb;
	int n1, n2;
	int iinfo = 0;
	/*
	*     Test the input parameters
	*/
	info = 0;
	bool upper = cblas_lsame(uplo, 'U');

	if (!upper & !cblas_lsame(uplo, 'L')) {
		info = -1;
	}
	else if (n < 0) {
		info = -2;
	}
	else if (lda < std::max(1, n)) {
		info = -4;
	}
	if (info != 0) {
		cblas_xerbla("DPOTRF2", -info);
		return;
	}

	MData a(n, n, A, lda);
	//a.show();

	/*
	*     Quick return if possible
	*/
	if (n == 0) {
		return;
	}
	/*
	* N = 1 case
	*/
	if (n == 1) {
		/*
		* Test for non - positive - definiteness
		*/
		if ((a[0][0] <= zero) || std::isnan(a[0][0])) {
			info = 1;
			return;
		}
		/*
		* Factor
		*/
		a[0][0] = sqrt(a[0][0]);
	}
	else {
		/*
		*Use recursive code
		*/
		n1 = n / 2;
		n2 = n - n1;
		/*
		*Factor A11
		*/
		clapack_dpotrf2(uplo, n1, A, lda, iinfo);
		if (iinfo != 0) {
			info = iinfo;
			return;
		}

		//
		//a.show();
		//


		/*
		*Compute the Cholesky factorization A = U * *T * U
		*/
		if (upper) {
			/*
			* Updateand scale A12
			*/
			cblas_dtrsm('L', 'U', 'T', 'N', n1, n2, one, A, lda, a.sub(0, n1), lda);

			//
			//a.show();
			//



			/*
			*Update and factor A22
			*/
			cblas_dsyrk(uplo, 'T', n2, n1, -one, a.sub(0, n1), lda, one, a.sub(n1, n1), lda);

			//
			//a.show();
			//


			clapack_dpotrf2(uplo, n2, a.sub(n1, n1), lda, iinfo);

			if (iinfo != 0) {
				info = iinfo + n1;
				return;
			}
		}
		else {
			/*
			*Compute the Cholesky factorization A = L * L * *T
			*/

			/*
			*Update and scale A21
			*/
			cblas_dtrsm('R', 'L', 'T', 'N', n2, n1, one, A, lda, a.sub(n1, 0), lda);

			/*
			*Update and factor A22
			*/
			cblas_dsyrk(uplo, 'N', n2, n1, -one, a.sub(n1, 0), lda, one, a.sub(n1, n1), lda);
			clapack_dpotrf2(uplo, n2, a.sub(n1, n1), lda, iinfo);

			if (iinfo != 0) {
				info = iinfo + n1;
				return;
			}
		}
	}
}


void clapack_dpotrf(char uplo, int n, double* A, int lda, int info)
{
	double ONE = 1.0;
	double ZERO = 0.0;

	int i, j, nb;
	int jb;

	j = 0;
	/*
	*     Test the input parameters
	*/
	info = 0;
	bool upper = cblas_lsame(uplo, 'U');

	if (!upper & !cblas_lsame(uplo, 'L')) {
		info = -1;
	}
	else if (n < 0) {
		info = -2;
	}
	else if (lda < std::max(1, n)) {
		info = -4;
	}
	if (info != 0) {
		cblas_xerbla("DPOTRF", -info);
		return;
	}

	/*
	*     Quick return if possible
	*/
	if (n == 0) {
		return;
	}

	MData a(n, n, A, lda);
	//a.setData(A);

	/*
	*     Determine the block size for this environment.
	*/
	nb = ilaenv(1, "DPOTRF", uplo, n, -1, -1, -1);
	if ((nb <= 1) || (nb >= n)) {
		/*
		* Use unblocked code.
		*/
		clapack_dpotrf2(uplo, n, A, lda, info);
	}
	else {
		/*
		* Use blocked code.
		*/
		if (upper) {
			/*
			* Compute the Cholesky factorization A = U**T*U.
			*/
			for (j = 0; j < n; j += nb) {
				/*
				*  Update and factorize the current diagonal block and test
				*  for non - positive - definiteness.
				*/
				jb = std::min(nb, n - j + 1);
				cblas_dsyrk('U', 'T', jb, j - 1, -ONE, &A[1 * lda + j], lda, ONE, &A[j * lda + j], lda);
				clapack_dpotrf2('U', jb, &A[j * lda + j], lda, info);
				if (info != 0) {
					goto MARK30;
				}
				if ((j + jb) <= n) {
					/*
					*  Compute the current block row.
					*/
					cblas_dgemm('T', 'N', jb, n - j - jb + 1, j - 1, -ONE, &A[1 * lda + j], lda, &A[1 * lda + j + jb], lda, ONE, &A[j * lda + j + jb], lda);
					cblas_dtrsm('L', 'U', 'T', 'N', jb, n - j - jb + 1, ONE, &A[j * lda + j], lda, &A[j * lda + j + jb], lda);
				}
			}
		}
		else {

			/*
			* Compute the Cholesky factorization A = L*L**T.
			*/
			for (j = 0; j < n; j += nb) {
				/*
				* Update and factorize the current diagonal block and test
				* for non - positive - definiteness.
				*/
				jb = std::min(nb, n - j + 1);
				cblas_dsyrk('L', 'N', jb, j - 1, -ONE, &A[j * lda + 1], lda, ONE, &A[j * lda + j], lda);
				clapack_dpotrf2('L', jb, &A[j * lda + j], lda, info);
				if (info != 0) {
					goto MARK30;
				}
				if ((j + jb) <= n) {
					/*
					* Compute the current block column.
					*/
					cblas_dgemm('N', 'T', n - j - jb + 1, jb, j - 1, -ONE, &A[(j + jb) * lda + 1], lda, &A[j * lda + 1], lda, ONE, &A[(j + jb) * lda + j], lda);
					cblas_dtrsm('R', 'L', 'T', 'N', n - j - jb + 1, jb, ONE, &A[j * lda + j], lda, &A[(j + jb) * lda + j], lda);
				}
			} 
		}
	}

MARK30:
	info = info + j - 1;

MARK40:
	return;
}


//
//
//{
//	double ONE = 1.0;
//	double ZERO = 0.0;
//
//	int i, j, nb;
//	int jb;
//	/*
//	*     Test the input parameters
//	*/
//	info = 0;
//	bool upper = cblas_lsame(uplo, 'U');
//	if (!upper &  !cblas_lsame(uplo, 'L')) {
//		info = -1;
//	}
//	else if (n < 0) {
//		info = -2;
//	}
//	else if (lda < std::max(1, n)) {
//		info = -4;
//	}
//	if (info != 0) {
//		xerbla("DPOTRF2", -info);
//		return;
//	}
//
//	/*
//	*     Quick return if possible
//	*/
//	if (n == 0) {
//		return;
//	}
//
//	/*
//	*     Determine the block size for this environment.
//	*/
//	nb = ilaenv(1, "DPOTRF", uplo, n, -1, -1, -1);
//
//
//	//  DEBUG -- BEGIN --
//	std::cout << "DPOTRF2: nb = " << nb << std::endl;
//	//  DEBUG -- END --
//
//
//	if ((nb <= 1) || (nb >= n)) {
//		/*
//		* Use unblocked code.
//		*/
//		cblas_dpotrf2(uplo, n, A, lda, info);
//	}
//	else {
//		/*
//		* Use blocked code.
//		*/
//		if (upper) {
//			/*
//			* Compute the Cholesky factorization A = U**T*U.
//			*/
//			for (j = 0; j < n; j += nb) {
//				/*
//				*  Update and factorize the current diagonal block and test
//				*  for non - positive - definiteness.
//				*/
//				jb = std::min(nb, n - j + 1);
//				cblas_dsyrk('U', 'T', jb, j - 1, -ONE, &A[1 * lda + j], lda, ONE, &A[j * lda + j], lda);
//				cblas_dpotrf2('U', jb, &A[j * lda + j], lda, info);
//				if (info != 0) {
//					goto MARK30;
//				}
//				if ((j + jb) <= n) {
//					/*
//					*  Compute the current block row.
//					*/
//					cblas_dgemm('T', 'N', jb, n - j - jb + 1, j - 1, -ONE, &A[1 * lda + j], lda, &A[1 * lda + j + jb], lda, ONE, &A[j * lda + j + jb], lda);
//					cblas_dtrsm('L', 'U', 'T', 'N', jb, n - j - jb + 1, ONE, &A[j * lda + j], lda, &A[j * lda + j + jb], lda);
//				}
//			}
//		}
//		else {
//			/*
//			* Compute the Cholesky factorization A = L*L**T.
//			*/
//			for (j = 0; j < n; j += nb) {
//				/*
//				* Update and factorize the current diagonal block and test
//				* for non - positive - definiteness.
//				*/
//				jb = std::min(nb, n - j + 1);
//				cblas_dsyrk('L', 'N', jb, j - 1, -ONE, &A[j * lda + 1], lda, ONE, &A[j * lda + j], lda);
//				cblas_dpotrf2('L', jb, &A[j * lda + j], lda, info);
//				if (info != 0) {
//					goto MARK30;
//				}
//				if ((j + jb) <= n) {
//					/*
//					* Compute the current block column.
//					*/
//					cblas_dgemm('N', 'T', n - j - jb + 1, jb, j - 1, -ONE, &A[(j + jb)*lda + 1], lda, &A[j * lda + 1], lda, ONE, &A[(j + jb) * lda + j], lda);
//					cblas_dtrsm('R', 'L', 'T', 'N', n - j - jb + 1, jb, ONE, &A[j * lda + j], lda, &A[(j + jb)*lda + j], lda);
//				}
//			}
//		}
//	}
//	goto MARK40;
//MARK30:
//	info = info + j - 1;
//
//MARK40:
//	return;
//}
//
//

/****
\brief  DPOTRS solves a system of linear equations A*X = B with a symmetric
 positive definite matrix A using the Cholesky factorization
 A = U**T*U or A = L*L**T computed by DPOTRF.
 @param[in]  UPLO
*/
void clapack_dpotrs(char UPLO, int n, int nrhs, double* A, int lda, double* B, int ldb, int INFO)
{
	INFO = 0;
	bool upper = cblas_lsame(UPLO, 'U');
	if (!upper && !cblas_lsame(UPLO, 'L')) {
		INFO = -1;
	}
	else  if (n < 0) {
		INFO = -2;
	}
	else  if (nrhs < 0) {
		INFO = -3;
	}
	else  if (lda < std::max(1, n)) {
		INFO = -5;
	}
	else  if (ldb < std::min(1, n)) {
		INFO = -7;
	}

	if (INFO != 0) {
		cblas_xerbla("DPOTRS", -INFO);
		return;
	}

	// Quick return if possible 
	if (n == 0 || nrhs == 0)
	{
		return;
	}

	double one = 1.0;

	if (upper) {
		/*
		  *        Solve A*X = B where A = U**T *U.
		 *
		*        Solve U**T *X = B, overwriting B with X.
		*/
		cblas_dtrsm('L', 'U', 'T', 'N', n, nrhs,
			one, A, lda, B, ldb);
		/*
		  *        Solve U*X = B, overwriting B with X.
		 */
		cblas_dtrsm('L', 'U', 'N', 'N', n,
			nrhs, one, A, lda, B, ldb);
	}
	else {
		/*
		  *        Solve A*X = B where A = L*L**T.
		 *
		  *        Solve L*X = B, overwriting B with X.
		 */
		cblas_dtrsm('L', 'L', 'N', 'N', n,
			nrhs, one, A, lda, B, ldb);
		/*
		 *        Solve L**T *X = B, overwriting B with X.
		 */
		cblas_dtrsm('L', 'L', 'T', 'N', n, nrhs,
			one, A, lda, B, ldb);
	}
	return;
}


/***
 DPOSV computes the solution to a real system of linear equations
	A * X = B,
 where A is an N-by-N symmetric positive definite matrix and X and B
 are N-by-NRHS matrices.

 The Cholesky decomposition is used to factor A as
	A = U**T* U,  if UPLO = 'U', or
	A = L * L**T,  if UPLO = 'L',
 where U is an upper triangular matrix and L is a lower triangular
 matrix.  The factored form of A is then used to solve the system of
 equations A * X = B.
 ***/
void clapack_dposv(char uplo, int n, int nrhs, double* a, int lda, double* b, int ldb, int& info)
{
	/*
	* Test the input parameters.
	*/
	info = 0;
	if (!cblas_lsame(uplo, 'U') && !cblas_lsame(uplo, 'L')) {
		info = -1;
	}
	else if (n < 0) {
		info = -2;
	}
	else if (nrhs < 0) {
		info = -3;
	}
	else if (lda < std::max(1, n)) {
		info = -5;
	}
	else if (ldb < std::min(1, n)) {
		info = -7;
	}
	if (info != 0) {
		cblas_xerbla("DPOSV ", -info);
		return;
	}
	/*
	* Compute the Cholesky factorization A = U**T*U or A = L*L**T.
	*/
	clapack_dpotrf(uplo, n, a, lda, info);
	if (info == 0) {
		/*
		* Solve the system A*X = B, overwriting B with X.
		*/
		clapack_dpotrs(uplo, n, nrhs, a, lda, b, ldb, info);
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

	for (int i = 0; i < k; i++) {
		/*
		*Generate elementary reflector H(i) to annihilate A(i + 1:m, i)
		*/
		//int sub_m = (std::min)(i, m - 1);
		//clapack_dlarfg(m - i, A[i][i], A.sub(i + 1, sub_m), lda, tau[i]); 
		int sub_m = (std::min)(i + 1, m - 1);
		clapack_dlarfg(m - i, A[i][i], A.sub(sub_m, i), lda, tau[i]);


		//dlarfg(m, A[0], A + incA, incA, tau);

		//
		//A.show();
		//

		if (i < (n - 1)) {
			/*
			*Apply H(i) to A(i:m, i + 1 : n) from the left
			*/
			clapack_dlarf1f('L', m - i, n - i - 1, A.sub(i, i), lda, tau[i], A.sub(i, i + 1), lda, work);


			// dlarf1f('L', m - i + 1, n - i, A, n, tau, A + 1, incA, work);

			//
			//A.show();
			//
		}
	}
	return;
}




/***
\brief DLACPY copies all or part of one two-dimensional array to another.
*/
void dlacpy(char uplo, int m, int n,
	double* pA, int lda, double* pB, int ldb)
{
	MData a(m, n, pA, ldb);
	MData b(m, n, pB, ldb);


	if (cblas_lsame(uplo, 'U'))
	{
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < std::min(j, m); i++)
			{
				b(i, j) = a(i, j);
			}
		}
	}
	else if (cblas_lsame(uplo, 'L'))
	{
		for (int j = 0; j < n; j++)
		{
			for (int i = j; i < m; i++)
			{
				b(i, j) = a(i, j);
			}
		}
	}
	else {
		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < m; i++)
			{
				b(i, j) = a(i, j);
			}
		}
	}
	return;
}

/***
 DLARFT forms the triangular factor T of a real block reflector H
 of order n, which is defined as a product of k elementary reflectors.

 If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;

 If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.

 If STOREV = 'C', the vector which defines the elementary reflector
 H(i) is stored in the i-th column of the array V, and

	H  =  I - V * T * V**T

 If STOREV = 'R', the vector which defines the elementary reflector
 H(i) is stored in the i-th row of the array V, and

	H  =  I - V**T * T * V

  The shape of the matrix V and the storage of the vectors which define
  the H(i) is best illustrated by the following example with n = 5 and
  k = 3. The elements equal to 1 are not stored.

  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':

			   V = (  1       )                 V = (  1 v1 v1 v1 v1 )
				   ( v1  1    )                     (     1 v2 v2 v2 )
				   ( v1 v2  1 )                     (        1 v3 v3 )
				   ( v1 v2 v3 )
				   ( v1 v2 v3 )

  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':

			   V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
				   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
				   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
				   (     1 v3 )
				   (        1 )


 @param[in] direct Specifies the order in which the elementary reflectors are
		   multiplied to form the block reflector:
		   = 'F': H = H(1) H(2) . . . H(k) (Forward)
		   = 'B': H = H(k) . . . H(2) H(1) (Backward)
 @param[in] storev  Specifies how the vectors which define the elementary
		   reflectors are stored (see also Further Details):
		   = 'C': columnwise
		   = 'R': rowwise
 @param[in] n The order of the block reflector H. N >= 0.
 @param[in] k The order of the triangular factor T (= the number of
			  elementary reflectors). K >= 1.
 @param[in] pV The matrix V. dimension
							   (LDV,K) if STOREV = 'C'
							   (LDV,N) if STOREV = 'R'
 @param[in] ldv  The leading dimension of the array V.
				 If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
 @param[in] tau  TAU is DOUBLE PRECISION array, dimension (K)
				 TAU(i) must contain the scalar factor of the elementary
				 reflector H(i).
 @param[in] pT   The k by k triangular factor T of the block reflector.
				 If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
				 lower triangular. The rest of the array is not used.
 @param[in] ldt The leading dimension of the array T. LDT >= K.
*/
void clapack_dlarft(char direct, char storev, int n, int k, double* pV, int ldv,
	double* tau, double* pT, int ldt)
{
	double one = 1.0;
	double zero = 0.0;
	int lastv = 0;
	int prevlastv = 0;
	int i, j;

	// Quick Return 
	if (n == 0)
	{
		return;
	}

	MData v(n, ldv, pV, ldv);
	MData t(n, ldt, pT, ldt);

	if (cblas_lsame(direct, 'F'))
	{
		// Forward 
		prevlastv = n;

		for (i = 0; i < k; i++) // DO i = 1, k
		{
			prevlastv = std::max(i, prevlastv);
			if (tau[i] == zero) {
				/*
				*H(i) = I
				*/
				for (j = 0; j < i; j++) // DO j = 1, i
				{
					t(j, i) = zero;
				}
			}
			else {
				/*
				*general case
				*/
				if (cblas_lsame(storev, 'C')) {
					// Skip any trailing zeros.
					for (lastv = n; lastv < i + 1; lastv += -1) {
						if (v(lastv, i) != zero)
						{
							return;
						}
					}
					for (j = 0; j < i; j++) // DO j = 1, i - 1
					{
						t(j, i) = -tau[i] * v(i, j);
					}
					j = std::min(lastv, prevlastv);


					//std::cout << " T: " << std::endl;
					//showMatrix_d(t.sub(0, i), i + 1, 1, ldt);


					/*
					*T(1:i - 1, i) := -tau(i) * V(i : j, 1 : i - 1)^T * V(i : j, i)
					*/
					//cblas_dgemv('T', j - i, i - 1, -tau[i],
					//	v.sub(i + 1, 1), ldv, v.sub(i + 1, i), 1, one,
					//	t.sub(1, i), 1);

					cblas_dgemv('T', j - i, i - 1, -tau[i],
						v.sub(i, 0), ldv, v.sub(i, i - 1), 1, one,
						t.sub(0, i - 1), 1);

					//cblas_dgemv('T', i + 1, j , -tau[i],
					//	v.sub( 0, 0), ldv, v.sub(0, i), ldv, one,
					//	t.sub(0, i), ldt);

					//cblas_dgemv('N', i - 1, j - i, -tau[i],
					//	v.sub(1, i + 1), ldv, v.sub(i, i + 1), ldv, one,
					//	t.sub(1, i), 1);

					//cblas_dgemv('T', i, j, -tauVec[i],
					//	A + i * lda, lda, A + i + j * lda, lda, one,
					//	T + i, lda);


					//std::cout << std::endl;
					//std::cout << " M: " << std::endl;
					//showMatrix_d(v.sub(0, 0), i + 1, i, ldv);
					//std::cout << " V: " << std::endl;
					//showMatrix_d(v.sub(0, i), i + 1, 1, ldv);
					//double tmp[]{ 0.0, 0.0, 0.0 };
					//cblas_dgemv('T', i + 1, i, -tau[i],
					//	v.sub(0, 0), ldv, v.sub(0, i), ldv, one,
					//	tmp, 1);
					//std::cout << " tmp: " << std::endl;
					//showMatrix_d(tmp, i + 1, 1, 1);


					//std::cout << " after tau(i) *  " << std::endl;
					//std::cout << " j: " << j << std::endl;
					//std::cout << " i: " << i << std::endl;
					//std::cout << " -tau[i]: " << -tau[i] << std::endl;
					//std::cout << " T: " << std::endl;
					//showMatrix_d(pT, 3, 3);

				}
				else {
					// Skip any trailing zeros.
					for (lastv = n; lastv < i + 1; lastv += -1) {
						if (v(lastv, i) != zero)
						{
							return;
						}
					}
					//DO lastv = n, i + 1, -1
					//	IF(v(i, lastv).NE.zero) EXIT
					//	END DO
					for (j = 0; j < i; j++) // DO j = 1, i - 1
					{
						t(j, i) = -tau[i] * v(j, i);
					}
					j = std::min(lastv, prevlastv);
					/*
					*T(1:i - 1, i) := -tau(i) * V(1:i - 1, i : j) * V(i, i:j)^T
					*/
					cblas_dgemv('N', i - 1, j - i, -tau[i],
						v.sub(0, i), ldv, v.sub(i - 1, i), ldv, one,
						t.sub(0, i - 1), 1);
				}
				/*
				*T(1:i - 1, i) := T(1:i - 1, 1 : i - 1) * T(1:i - 1, i)
				*/
				cblas_dtrmv('U', 'N', 'N', i, t._pData, ldt, t.sub(0, i), ldt);
				t(i, i) = tau[i];

				////
				//std::cout << " after trmv  T:" << std::endl;
				//showMatrix_d(pT, 3, 3);
				////

				if (i > 1) {
					prevlastv = std::max(prevlastv, lastv);
				}
				else {
					prevlastv = lastv;
				}
			}
		} //  
	}
	else {
		// Backward 
		prevlastv = 0;
		for (i = k - 1; i >= 0; i--)
		{
			if (tau[i] == zero)
			{
				/*
				*H(i) = I
				*/
				for (j = i - 1; j < k; j++) // DO j = i, k
				{
					t(j, i) = zero;
				}
			}
			else {
				/*
				*general case
				*/
				if (cblas_lsame(storev, 'C'))
				{
					// Skip any trailing zeros. 
					for (lastv = 0; lastv < i - 1; lastv++) {
						if (v(lastv, i) != zero)
						{
							return;
						}
					}
					for (j = i; j < k; j++) // j = i+1, k
					{
						t(j, i) = -tau[i] * v(n - k + i, j);
					}
					j = std::max(lastv, prevlastv);
					/*
					 * T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**T * V(j:n-k+i,i)
					 */
					cblas_dgemv('T', n - k + i - j, k - i, -tau[i],
						v.sub(j - 1, i), ldv, v.sub(j - 1, i - 1), 1, one,
						t.sub(i, i - 1), 1);
				}
				else {
					// Skip any leading zeros.
					for (lastv = 0; lastv < i; lastv++)    //DO lastv = 1, i - 1
					{
						if (v(i, lastv) != zero)
						{
							return;
						}
					} // END DO
					for (j = i; j < k; j++) // //  DO j = i + 1, k
					{
						t(j, i) = -tau[i] * v(j, n - k + i);
					} // END DO 
					j = std::max(lastv, prevlastv);
					/*
					*T(i + 1:k, i) = -tau(i) * V(i + 1:k, j : n - k + i) * V(i, j:n - k + i) * *T
					*/
					cblas_dgemv('N', k - i, n - k + i - j,
						-tau[i], v.sub(i, j - 1), ldv, v.sub(i - 1, j - 1), ldv,
						one, t.sub(i, i - 1), 1);
				} // EMD IF  
				  /*
				  *  T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
				  */
				cblas_dtrmv('L', 'N', 'N', k - i, t.sub(i + 1, i + 1), ldt, t.sub(i + 1, i), 1);
				if (i > 0) {
					prevlastv = std::min(prevlastv, lastv);
				}
				else {
					prevlastv = lastv;
				}
			}
			t(i, i) = tau[i];
		}
	}
}


/***
 DLARFB applies a real block reflector H or its transpose H^T to a
 real m by n matrix C, from either the left or the right.

  The shape of the matrix V and the storage of the vectors which define
  the H(i) is best illustrated by the following example with n = 5 and
  k = 3. The triangular part of V (including its diagonal) is not
  referenced.

  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':

			   V = (  1       )                 V = (  1 v1 v1 v1 v1 )
				   ( v1  1    )                     (     1 v2 v2 v2 )
				   ( v1 v2  1 )                     (        1 v3 v3 )
				   ( v1 v2 v3 )
				   ( v1 v2 v3 )

  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':

			   V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
				   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
				   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
				   (     1 v3 )
				   (        1 )

*/
void clapack_dlarfb(char side, char trans, char direct, char storev,
	int m, int n, int k, double* v, int ldv, double* t,
	int ldt, double* pC, int ldc, double* work, int ldwork)
{
	/*
	*     Quick return; if possible
	*/
	if ((m <= 0) || (n <= 0))
	{
		return;
	}

	int i, j;
	char transt;
	double one = 1.0;

	MData V(m, n, v, ldv);
	MData c(m, n, pC, ldc);
	MData Work(m, n, work, ldwork);

	if (cblas_lsame(trans, 'N')) {
		transt = 'T';
	}
	else {
		transt = 'N';
	}



	if (cblas_lsame(storev, 'C')) {

		if (cblas_lsame(direct, 'F')) {
			/*
			*           Let  V =  ( V1 )    (first K rows)
			*                     ( V2 )
			*           where  V1  is U L triangular.
			*/
			if (cblas_lsame(side, 'L')) {
				/*
				*              Form  H * C  or  H**T * C  where  C = ( C1 )
				*                                                    ( C2 )
				*
				*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
				*
				*              W := C1**T
				*/
				for (j = 0; j < k; j++) {
					cblas_dcopy(n, c.sub(0, j), ldc, Work.sub(j, 0), 1);
				}

				//
				//printf("   dlarfb() :  W := C1**T  \n");
				//Work.show();
				//

				/*
				*              W := W * V1 = C^T V
				*/
				cblas_dtrmm('R', 'L', 'N', 'U',
					n, k, one, v, ldv, work, ldwork);


				//
				//printf("   dlarfb() :  W := W * V1 = C^T V  \n");
				//Work.show();
				//

				if (m > k) {
					/*
					*                 W := W + C2**T * V2
					*/
					cblas_dgemm('T', 'N', n, k, m - k,
						one, c.sub(k, 0), ldc, V.sub(k, 0), ldv,
						one, work, ldwork);
				}
				/*
				*              W := W * T**T  or  W * T
				*/
				cblas_dtrmm('R', 'U', transt, 'N', n,
					k, one, t, ldt, work, ldwork);
				/*
				*              C := C - V * W**T
				*/
				if (m > k) {
					/*
					*                 C2 := C2 - V2 * W**T
					*/
					cblas_dgemm('N', 'T', m - k, n, k,
						-one, V.sub(k, 0), ldv, work, ldwork, one,
						c.sub(k, 0), ldc);
				}
				/*
				*              W := W * V1**T
				*/
				cblas_dtrmm('R', 'L', 'T', 'U', n,
					k, one, v, ldv, work, ldwork);
				/*
				*              C1 := C1 - W**T
				*/
				for (j = 0; j < k; j++) {
					for (i = 0; i < n; i++) {
						c(j, i) = c(j, i) - Work(i, j);
					}
				}

			}
			else if (cblas_lsame(side, 'R')) {
				/*
				*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
				*
				*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
				*
				*              W := C1
				*/
				for (j = 0; j < k; j++) {
					cblas_dcopy(m, c.sub(0, j), 1, Work.sub(0, j), 1);
				}
				/*
				*              W := W * V1
				*/
				cblas_dtrmm('R', 'L', 'N', 'U',
					m, k, one, v, ldv, work, ldwork);
				if (n > k) {
					/*
					*                 W := W + C2 * V2
					*/
					cblas_dgemm('N', 'N', m, k,
						n - k, one, c.sub(0, k), ldc, V.sub(k, 0), ldv,
						one, work, ldwork);
				}
				/*
				*              W := W * T  or  W * T**T
				*/
				cblas_dtrmm('R', 'U', trans, 'N', m, k,
					one, t, ldt, work, ldwork);
				/*
				*   C := C - W * V**T = C - C V T V^T
				*     = C(I - VTV^T)
				*/
				if (n > k) {
					/*
					*                 C2 := C2 - W * V2**T
					*/
					cblas_dgemm('N', 'T', m, n - k, k,
						-one, work, ldwork, V.sub(k, 0), ldv, one,
						c.sub(0, k), ldc);
				}
				/*
				*              W := W * V1**T
				*/
				cblas_dtrmm('R', 'L', 'T', 'U', m,
					k, one, v, ldv, work, ldwork);
				/*
				*              C1 := C1 - W
				*/
				for (j = 0; j < k; j++) {
					for (i = 0; i < m; i++) {
						c(i, j) = c(i, j) - Work(i, j);
					}
				}
			}

		}
		else {
			/*
			*           Let  V =  ( V1 )
			*                     ( V2 )    (last K rows)
			*           where  V2  is U U triangular.
			*/
			if (cblas_lsame(side, 'L')) {
				/*
				*              Form  H * C  or  H**T * C  where  C = ( C1 )
				*                                                    ( C2 )
				*
				*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
				*
				*              W := C2**T
				*/
				for (j = 0; j < k; j++) {
					cblas_dcopy(n, c.sub(0, m - k + j), ldc, Work.sub(j, 0), 1);
				}

				//
				printf("   dlarfb()  W = C2^T   \n");
				Work.show();
				//

				/*
				*              W := W * V2 = C2^T * V2
				*/
				cblas_dtrmm('R', 'U', 'N', 'U',
					n, k, one, V.sub(m - k, 0), ldv, work, ldwork);

				//
				printf("   dlarfb()   W := W * V2 = C2^T * V2   \n");
				Work.show();
				//

				if (m > k) {
					/*
					*     W := W + C1**T * V1
					*/
					cblas_dgemm('T', 'N', n, k, m - k,
						one, pC, ldc, v, ldv, one, work, ldwork);
				}
				/*
				*      W := W * T**T  or  W * T
				*/
				cblas_dtrmm('R', 'L', transt, 'N', n,
					k, one, t, ldt, work, ldwork);

				//
				printf("   dlarfb()   W := W * T**T    \n");
				Work.show();
				//

				/*
				*      C := C - V * W**T
				*/
				if (m > k) {
					/*
					*                 C1 := C1 - V1 * W**T
					*/
					cblas_dgemm('N', 'T', m - k, n, k,
						-one, v, ldv, work, ldwork, one, pC, ldc);
				}
				/*
				*              W := W * V2**T
				*/
				cblas_dtrmm('R', 'U', 'T', 'U', n,
					k, one, V.sub(m - k + 1, 1), ldv, work, ldwork);
				/*
				*              C2 := C2 - W**T
				*/
				for (j = 0; j < k; j++) {
					for (i = 0; i < n; i++) {
						c(m - k + j, i) = c(m - k + j, i) - Work(i, j);
					}
				}

			}
			else if (cblas_lsame(side, 'R')) {
					/*
					*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
					*
					*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
					*
					*              W := C2
					*/
					for (j = 0; j < k; j++) {
						cblas_dcopy(m, c.sub(0, n - k + j), 1, Work.sub(0, j), 1);
					}

					//
					Work.show();
					//

					/*
					*              W := W * V2
					*/
					cblas_dtrmm('R', 'U', 'N', 'U',
						m, k, one, V.sub(n - k, 0), ldv, work, ldwork);
					if (n > k) {
						/*
						*                 W := W + C1 * V1
						*/
						cblas_dgemm('N', 'N', m, k,
							n - k, one, pC, ldc, v, ldv, one, work, ldwork);
					}
					/*
					*              W := W * T  or  W * T**T
					*/
					cblas_dtrmm('R', 'L', trans, 'N', m, k,
						one, t, ldt, work, ldwork);
					/*
					*              C := C - W * V**T
					*/
					if (n > k) {
						/*
						*                 C1 := C1 - W * V1**T
						*/
						cblas_dgemm('N', 'T', m, n - k, k,
							-one, work, ldwork, v, ldv, one, pC, ldc);
					}
					/*
					*              W := W * V2**T
					*/
					cblas_dtrmm('R', 'U', 'T', 'U', m,
						k, one, V.sub(n - k, 0), ldv, work, ldwork);
					/*
					*              C2 := C2 - W
					*/
					for (j = 0; j < k; j++) {
						for (i = 0; i < m; i++) {
							c(i, n - k + j) = c(i, n - k + j) - Work(i, j);
						}
					}
				}
		}

	}
	else if (cblas_lsame(storev, 'R')) {

		if (cblas_lsame(direct, 'F')) {
			/*
			*           Let  V =  ( V1  V2 )    (V1: first K columns)
			*           where  V1  is U U triangular.
			*/
			if (cblas_lsame(side, 'L')) {
				/*
				*              Form  H * C  or  H**T * C  where  C = ( C1 )
				*                                                    ( C2 )
				*
				*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
				*
				*              W := C1**T
				*/
				for (j = 0; j < k; j++) {
					cblas_dcopy(n, c.sub(0, j), ldc, Work.sub(j, 0), 1);
				}
				//
				printf("   dlarfb()  W = C1^T   \n");
				Work.show();
				//

				/*
				*              W := W * V1**T
				*/
				cblas_dtrmm('R', 'U', 'T', 'U', n,
					k, one, v, ldv, work, ldwork);

				//
				printf("   dlarfb()  W := W * V1**T = C1^T * V1**T    \n");
				Work.show();
				printf("             V =    \n");
				V.show();
				//

				if (m > k) {
					/*
					*                 W := W + C2**T * V2**T
					*/
					cblas_dgemm('T', 'T', n, k, m - k,
						one, c.sub(k, 0), ldc, V.sub(0, k), ldv, one,
						work, ldwork);
				}
				/*
				*              W := W * T**T  or  W * T
				*/
				cblas_dtrmm('R', 'U', transt, 'N', n,
					k, one, t, ldt, work, ldwork);
				/*
				*              C := C - V**T * W**T
				*/
				if (m > k) {
					/*
					*                 C2 := C2 - V2**T * W**T
					*/
					cblas_dgemm('T', 'T', m - k, n, k,
						-one, V.sub(0, k), ldv, work, ldwork, one,
						c.sub(k, 0), ldc);
				}
				/*
				*              W := W * V1
				*/
				cblas_dtrmm('R', 'U', 'N', 'U',
					n, k, one, v, ldv, work, ldwork);
				/*
				*              C1 := C1 - W**T
				*/
				for (j = 0; j < k; j++) {
					for (i = 0; i < n; j++) {
						c(j, i) = c(j, i) - Work(i, j);
					}
				}

			}
			else if (cblas_lsame(side, 'R')) {
					/*
					*              Form  C * H  or  C * H**T  where  C = ( C1  C2 )
					*
					*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
					*
					*              W := C1
					*/
					for (j = 0; j < k; j++) {
						cblas_dcopy(m, c.sub(0, j), 1, Work.sub(0, j), 1);
					}
					/*
					*              W := W * V1**T
					*/
					cblas_dtrmm('R', 'U', 'T', 'U', m,
						k, one, v, ldv, work, ldwork);
					if (n > k) {
						/*
						*                 W := W + C2 * V2**T
						*/
						cblas_dgemm('N', 'T', m, k, n - k,
							one, c.sub(0, k), ldc, V.sub(0, k), ldv,
							one, work, ldwork);
					}
					/*
					*              W := W * T  or  W * T**T
					*/
					cblas_dtrmm('R', 'U', trans, 'N', m, k,
						one, t, ldt, work, ldwork);
					/*
					*              C := C - W * V
					*/
					if (n > k) {
						/*
						*                 C2 := C2 - W * V2
						*/
						cblas_dgemm('N', 'N', m, n - k,
							k, -one, work, ldwork, V.sub(0, k), ldv, one,
							c.sub(0, k), ldc);
					}
					/*
					*              W := W * V1
					*/
					cblas_dtrmm('R', 'U', 'N', 'U',
						m, k, one, v, ldv, work, ldwork);
					/*
					*              C1 := C1 - W
					*/
					for (j = 0; j < k; j++) {
						for (i = 0; i < m; i++) {
							c(i, j) = c(i, j) - Work(i, j);
						}
					}

				}

		}
		else {
			/*
			*           Let  V =  ( V1  V2 )    (V2: last K columns)
			*           where  V2  is U L triangular.
			*/
			if (cblas_lsame(side, 'L')) {
				/*
				*              Form  H * C  or  H**T * C  where  C = ( C1 )
				*                                                    ( C2 )
				*
				*              W := C**T * V**T  =  (C1**T * V1**T + C2**T * V2**T) (stored in WORK)
				*
				*              W := C2**T
				*/
				for (j = 0; j < k; j++) {
					cblas_dcopy(n, c.sub(0, m - k + j), ldc, Work.sub(j, 0), 1);
				}

				//
				printf("   dlarfb()  Work: W := C2**T  \n");
				Work.show();
				//

				/*
				*              W := W * V2**T
				*/
				cblas_dtrmm('R', 'L', 'T', 'U', n,
					k, one, V.sub(0, m - k), ldv, work, ldwork);
				if (m > k) {
					/*
					*                 W := W + C1**T * V1**T
					*/
					cblas_dgemm('T', 'T', n, k, m - k,
						one, pC, ldc, v, ldv, one, work, ldwork);
				}

				//
				printf("   dlarfb()  Work: W := W * V2**T = C2**T * V2**T \n");
				Work.show();
				//

				/*
				*              W := W * T**T  or  W * T
				*/
				cblas_dtrmm('R', 'L', transt, 'N', n,
					k, one, t, ldt, work, ldwork);
				/*
				*              C := C - V**T * W**T
				*/
				if (m > k) {
					/*
					*                 C1 := C1 - V1**T * W**T
					*/
					cblas_dgemm('T', 'T', m - k, n, k,
						-one, v, ldv, work, ldwork, one, pC, ldc);
				}
				/*
				*              W := W * V2
				*/
				cblas_dtrmm('R', 'L', 'N', 'U',
					n, k, one, V.sub(0, m - k), ldv, work, ldwork);
				/*
				*              C2 := C2 - W**T
				*/
				for (j = 0; j < k; j++) {
					for (i = 0; i < n; i++) {
						c(m - k + j, i) = c(m - k + j, i) - Work(i, j);
					}
				}

			}
			else if (cblas_lsame(side, 'R')) {
				/*
				*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
				*
				*              W := C * V**T  =  (C1*V1**T + C2*V2**T)  (stored in WORK)
				*
				*              W := C2
				*/
				for (j = 0; j < k; j++) {
					cblas_dcopy(m, c.sub(0, n - k + j), 1, Work.sub(0, j), 1);
				}
				/*
				*              W := W * V2**T
				*/
				cblas_dtrmm('R', 'L', 'T', 'U', m,
					k, one, V.sub(1, n - k), ldv, work, ldwork);
				if (n > k) {
					/*
					*                 W := W + C1 * V1**T
					*/
					cblas_dgemm('N', 'T', m, k, n - k,
						one, pC, ldc, v, ldv, one, work, ldwork);
				}
				/*
				*              W := W * T  or  W * T**T
				*/
				cblas_dtrmm('R', 'L', trans, 'N', m, k,
					one, t, ldt, work, ldwork);
				/*
				*              C := C - W * V
				*/
				if (n > k) {
					/*
					*                 C1 := C1 - W * V1
					*/
					cblas_dgemm('N', 'N', m, n - k,
						k, -one, work, ldwork, v, ldv, one, pC, ldc);
				}
				/*
				*              W := W * V2
				*/
				cblas_dtrmm('R', 'L', 'N', 'U',
					m, k, one, V.sub(1, n - k), ldv, work, ldwork);
				/*
				*              C1 := C1 - W
				*/
				for (j = 0; j < k; j++) {
					for (i = 0; i < m; i++) {
						c(i, n - k + j) = c(i, n - k + j) - Work(i, j);
					}
				}

			}
		}
	}
	return;
}


/***
DGEQRF computes a QR factorization of a real M-by-N matrix A:

	A = Q * ( R ),
			( 0 )

 where:

	Q is a M-by-M orthogonal matrix;
	R is an upper-triangular N-by-N matrix;
	0 is a (M-N)-by-N zero matrix, if M > N.

  The matrix Q is represented as a product of elementary reflectors

	 Q = H(1) H(2) . . . H(k), where k = min(m,n).

  Each H(i) has the form

	 H(i) = I - tau * v * v**T

  where tau is a real scalar, and v is a real vector with
  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
  and tau in TAU(i).

\param[in] m: int, The number of rows of the matrix A.  M >= 0.
\param[in] n: int, The number of columns of the matrix A.  N >= 0.
\param[in] a: double*, dimension (LDA,N)
		  On entry, the M-by-N matrix A.
		  On exit, the elements on and above the diagonal of the array
		  contain the min(M,N)-by-N upper trapezoidal matrix R (R is
		  upper triangular if m >= n); the elements below the diagonal,
		  with the array TAU, represent the orthogonal matrix Q as a
		  product of min(m,n) elementary reflectors (see Further
		  Details).
\param[in] lda: int, The leading dimension of the array A.  LDA >= max(1,M).
\param[in] tau: double*, dimension (min(M,N))
		  The scalar factors of the elementary reflectors (see Further
		  Details).
\param[in] work: double*, dimension (MAX(1,LWORK))
		  On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
\param[in] lowork: int, The dimension of the array WORK.
		  LWORK >= 1, if MIN(M,N) = 0, and LWORK >= N, otherwise.
		  For optimum performance LWORK >= N*NB, where NB is
		  the optimal blocksize.

		  If LWORK = -1, then a workspace query is assumed; the routine
		  only calculates the optimal size of the WORK array, returns
		  this value as the first entry of the WORK array, and no error
		  message related to LWORK is issued by XERBLA.
\param[out] info: int
		  = 0:  successful exit
		  < 0:  if INFO = -i, the i-th argument had an illegal value
*/
void clapack_dgeqrf(int	m, int n, double* pA, int lda, double* tau, double* work, int lwork, int& info)
{
	int iinfo;
	int lwkopt;
	int i;
	/*
	* Test the input arguments
	*/
	int k = (std::min)(m, n);
	info = 0;
	int nb = ilaenv(1, "DGEQRF", ' ', m, n, -1, -1);
	bool lquery = (lwork == -1);


	if (m < 0) {
		info = -1;
	}
	else if (n < 0) {
		info = -2;
	}
	else if (lda < (std::max)(1, n)) {
		info = -4;
	}
	else if (!lquery) {
		if ((lwork <= 0) || ((m > 0) && lwork < (std::max)(1, n))) {
			info = -7;
		}
	}
	if (info != 0) {
		cblas_xerbla("DGEQRF", -info);
		return;
	}
	else if (lquery) {
		//If LWORK = -1, then a workspace query is assumed; the routine
		//	only calculates the optimal size of the WORK array, returns
		//	this value as the first entry of the WORK array, and no error
		//	message related to LWORK is issued by XERBLA.

		if (k == 0) {
			lwkopt = 1;
		}
		else {
			lwkopt = n * nb;
		}
		work[0] = lwkopt;
		return;
	}

	MData a(m, n, pA, lda);

	/*
	*Quick return if possible
	*/
	if (k == 0) {
		work[0] = 0;
		return;
	}

	//int lwork;
	int ldwork;
	int nbmin = 2;
	int nx = 0;
	int iws = n;

	//printf("   nb: %i   \n", nb);
	//printf("    k: %i   \n", k);


	if (nb > 1 && nb < k) {
		/*
		*  Determine when to cross over from blocked to unblocked code.
		*/
		nx = (std::max)(0, ilaenv(3, "DGEQRF", ' ', m, n, -1, -1));
		if (nx < k) {
			/*
			*   Determine if workspace is large enough for blocked code.
			*/
			ldwork = n;
			iws = ldwork * nb;
			if (lwork < iws) {
				/*
				*  Not enough workspace to use optimal NB:  reduce NB and
				*  determine the minimum value of NB.
				*/
				nb = lwork / ldwork;
				nbmin = (std::max)(2, ilaenv(2, "DGEQRF", ' ', m, n, -1, -1));
			}
		}
	}
	if ((nb >= nbmin) && (nb < k) && (nx < k)) {
		/*
		*        Use blocked code initially
		*/
		for (i = 0; i < k - nx; i += nb) {
			int ib = (std::min)(k - i, nb);
			/*
			*           Compute the QR factorization of the current block
			*           A(i:m,i:i+ib-1)
			*/

			//printf("     \n");
			//printf("   Compute the QR factorization of the current block  \n");
			//printf("    A(i:m,i:i+ib-1)   \n");
			//printf("     \n");

			clapack_dgeqr2(m - i, ib, a.sub(i, i), lda, tau + i, work, iinfo);


			if (i + ib <= n) {
				/*
				*   Form the triangular factor of the block reflector
				*   H = H(i) H(i+1) . . . H(i+ib-1)
				*/
				clapack_dlarft('F', 'C', m - i, ib, a.sub(i, i), lda, tau + i, work, ldwork);
				/*
				*  Apply H**T to A(i:m,i+ib:n) from the left
				*/
				clapack_dlarfb('L', 'T', 'F', 'C', m - i, n - i - ib, ib,
					a.sub(i, i), lda, work, ldwork, a.sub(i, i + ib),
					lda, work + (ib), ldwork);
			}
		}
	}
	else {
		i = 0;
	}
	/*
	*     Use unblocked code to factor the last or only block.
	*/
	if (i <= k) {
		clapack_dgeqr2(m - i, n - i, a.sub(i, i), lda, tau + i, work, iinfo);
	}
	work[0] = iws;
	return;
}






/******
 DLARGV generates a vector of plane rotations with real cosines and real sines.
 DLARGV generates a vector of real plane rotations, determined by
 elements of the real vectors x and y. For i = 1,2,...,n

	(  c(i)  s(i) ) ( x(i) ) = ( a(i) )
	( -s(i)  c(i) ) ( y(i) ) = (   0  )
*****/
void dlargv(int n, double* x, int incx, double* y, int incy, double* c, int incc)
{
	double f, g, t, tt;
	int ix, iy, ic;
	double zero = 0.0;
	double one = 1.0;
	ix = 0;
	iy = 0;
	ic = 0;
	for (int i = 0; i < n; i++) {
		f = x[ix];
		g = y[iy];
		if (g == zero) {
			c[ic] = one;
		}
		else if (f == zero) {
			c[ic] = zero;
			y[iy] = one;
			x[ix] = g;
		}
		else if (fabs(f) > fabs(g)) {
			t = g / f;
			tt = sqrt(one + t * t);
			c[ic] = one / tt;
			y[iy] = t * c[ic];
			x[ix] = f * tt;
		}
		else {
			t = f / g;
			tt = sqrt(one + t * t);
			y[iy] = one / tt;
			c[ic] = t * y[iy];
			x[ix] = g * tt;
		}
		ic = ic + incc;
		iy = iy + incy;
		ix = ix + incx;
	}
	return;
}

/***
 DLARTV applies a vector of real plane rotations to elements of the
 real vectors x and y. For i = 1,2,...,n

	( x(i) ) := (  c(i)  s(i) ) ( x(i) )
	( y(i) )    ( -s(i)  c(i) ) ( y(i) )
 ****/
void dlartv(int n, double* x, int incx, double* y, int incy, double* c, double* s, int incc)
{
	double xi, yi;
	int ix = 0;
	int iy = 0;
	int ic = 0;
	for (int i = 1; i < n; i++) {
		xi = x[ix];
		yi = y[iy];
		x[ix] = c[ic] * xi + s[ic] * yi;
		y[iy] = c[ic] * yi - s[ic] * xi;
		ix = ix + incx;
		iy = iy + incy;
		ic = ic + incc;
	}
	return;
}