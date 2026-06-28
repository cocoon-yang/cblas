#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include "cblas.h"
#include <algorithm>
#include "mdata.h"

//#define min(x,y) (((x) < (y)) ? (x) : (y))
//#define max(x,y) (((x) > (y)) ? (x) : (y)) 

extern double cblas_idamax(int n, double* dx, int incx);
extern double dlamch(const char* CMACH);
extern int ilaenv(int ispec, const char* name, char OPTS, int N1, int N2, int N3, int N4);
extern bool cblas_lsame(char CA, char CB);
extern void cblas_xerbla(const char* SRNAME, int INFO);
extern void cblas_dscal(const int N, const double DA, double *DX, const int INCX);
extern void cblas_dlaswp(int n, double* pA, int lda, int k1, int k2, int* ipiv, int incx); 

//==============================================================
// Level 3

/******************************************************************************/

/******************************************************************************/
/*
  Purpose:

    DGEMM computes C = alpha * A * B and related operations.

  Discussion:

    DGEMM performs one of the matrix-matrix operations

     C := alpha * op ( A ) * op ( B ) + beta * C,

    where op ( X ) is one of

      op ( X ) = X   or   op ( X ) = X',

    ALPHA and BETA are scalars, and A, B and C are matrices, with op ( A )
    an M by K matrix, op ( B ) a K by N matrix and C an N by N matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.
    
  Modified:

    10 February 2014

  Author:

    Original FORTRAN77 version by Jack Dongarra.
    C version by John Burkardt.

  Parameters:

    Input, char TRANSA, specifies the form of op( A ) to be used in
    the matrix multiplication as follows:
    'N' or 'n', op ( A ) = A.
    'T' or 't', op ( A ) = A'.
    'C' or 'c', op ( A ) = A'.

    Input, char TRANSB, specifies the form of op ( B ) to be used in
    the matrix multiplication as follows:
    'N' or 'n', op ( B ) = B.
    'T' or 't', op ( B ) = B'.
    'C' or 'c', op ( B ) = B'.

    Input, int M, the number of rows of the  matrix op ( A ) and of the  
    matrix C.  0 <= M.

    Input, int N, the number  of columns of the matrix op ( B ) and the 
    number of columns of the matrix C.  0 <= N.

    Input, int K, the number of columns of the matrix op ( A ) and the 
    number of rows of the matrix op ( B ).  0 <= K.

    Input, double ALPHA, the scalar multiplier 
    for op ( A ) * op ( B ).

    Input, double A(LDA,KA), where:
    if TRANSA is 'N' or 'n', KA is equal to K, and the leading M by K
    part of the array contains A;
    if TRANSA is not 'N' or 'n', then KA is equal to M, and the leading
    K by M part of the array must contain the matrix A.

    Input, int LDA, the first dimension of A as declared in the calling 
    routine.  When TRANSA = 'N' or 'n' then LDA must be at least max ( 1, M ), 
    otherwise LDA must be at least max ( 1, K ).

    Input, double B(LDB,KB), where:
    if TRANSB is 'N' or 'n', kB is N, and the leading K by N 
    part of the array contains B;
    if TRANSB is not 'N' or 'n', then KB is equal to K, and the leading
    N by K part of the array must contain the matrix B.

    Input, int LDB, the first dimension of B as declared in the calling 
    routine.  When TRANSB = 'N' or 'n' then LDB must be at least max ( 1, K ), 
    otherwise LDB must be at least max ( 1, N ).

    Input, double BETA, the scalar multiplier for C.

    Input/output, double C(LDC,N).
    Before entry, the leading M by N part of this array must contain the 
    matrix C, except when BETA is 0.0, in which case C need not be set 
    on entry.
    On exit, the array C is overwritten by the M by N matrix
      alpha * op ( A ) * op ( B ) + beta * C.

    Input, int LDC, the first dimension of C as declared in the calling 
    routine.  max ( 1, M ) <= LDC.
*/
void cblas_dgemm(char transa, char transb, int m, int n, int k,
	double alpha, double* pA, int lda, double* pB, int ldb, double beta,
	double* pC, int ldc)
{
	int i;
	int j;
	int l;
	int nrowa;
	int nrowb;
	int ncola;
	int ncolb;

	int nota;
	int notb;
	double temp;
	/*
	Set NOTA and NOTB as true if A and B respectively are not
	transposed and set NROWA, NCOLA and NROWB as the number of rows
	and columns of A and the number of rows of B respectively.
	*/
	nota = ((transa == 'N') || (transa == 'n'));

	if (nota)
	{
		nrowa = m;
		ncola = k;
	}
	else
	{
		nrowa = k;
		ncola = m;
	}

	notb = ((transb == 'N') || (transb == 'n'));

	if (notb)
	{
		nrowb = k;
		ncolb = n;
	}
	else
	{
		nrowb = n;
		ncolb = k;
	}
	/*
	Test the input parameters.
	*/
	if (!(transa == 'N' || transa == 'n' ||
		transa == 'C' || transa == 'c' ||
		transa == 'T' || transa == 't'))
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "DGEMM - Fatal error!\n");
		fprintf(stderr, "  Input TRANSA had illegal value.\n");
		//exit ( 1 );
		return;
	}

	if (!(transb == 'N' || transb == 'n' ||
		transb == 'C' || transb == 'c' ||
		transb == 'T' || transb == 't'))
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "DGEMM - Fatal error!\n");
		fprintf(stderr, "  Input TRANSB had illegal value.\n");
		//exit ( 1 );
		return;
	}

	if (m < 0)
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "DGEMM - Fatal error!\n");
		fprintf(stderr, "  Input M had illegal value.\n");
		//exit ( 1 );
		return;
	}

	if (n < 0)
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "DGEMM - Fatal error!\n");
		fprintf(stderr, "  Input N had illegal value.\n");
		//exit ( 1 );
		return;
	}

	if (k < 0)
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "DGEMM - Fatal error!\n");
		fprintf(stderr, "  Input K had illegal value.\n");
		//exit ( 1 );
		return;
	}

	//if ( lda < i4_max ( 1, nrowa ) )
	if (lda < std::max(1, ncola))
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "DGEMM - Fatal error!\n");
		fprintf(stderr, "  Input LDA had illegal value.\n");
		//exit ( 1 );
		return;
	}

	if (ldb < std::max(1, ncolb))
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "DGEMM - Fatal error!\n");
		fprintf(stderr, "  Input LDB had illegal value.\n");
		//exit ( 1 );
		return;
	}

	if (ldc < std::max(1, n))//if ( ldc < i4_max ( 1, m ) )
	{
		fprintf(stderr, "\n");
		fprintf(stderr, "DGEMM - Fatal error!\n");
		fprintf(stderr, "  Input LDC had illegal value.\n");
		//exit ( 1 );
		return;
	}
	/*
	Quick return if possible.
	*/
	if (m == 0)
	{
		return;
	}

	if (n == 0)
	{
		return;
	}

	if ((alpha == 0.0 || k == 0) && (beta == 1.0))
	{
		return;
	}

	MData a(m, lda, pA, lda);
	//a.setData(pA);

	MData b(k, ldb, pB, ldb);
	//b.setData(pB);

	MData c(m, ldc, pC, ldc);
	//c.setData(pC);

	/*
	And if alpha is 0.0.
	*/
	if (alpha == 0) {
		if (beta == 0) {
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < m; i++) {
					c[i][j] = 0;
				}
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < m; i++) {
					c[i][j] = beta * c[i][j];
				}
			}
		}
		return;
	}

	/*
	Start the operations.
	*/
	if (notb)
	{
		if (nota)
		{
			for (int j = 0; j < n; j++) {
				if (beta == 0) {
					for (int i = 0; i < m; i++) {
						c(i, j) = 0;
					}
				}
				else if (beta != 1) {
					for (int i = 0; i < m; i++) {
						c(i, j) = beta * c(i, j);
					}
				}
				for (int l = 0; l < k; l++) {
					temp = alpha * b(l, j);
					for (int i = 0; i < m; i++) {
						c(i, j) = c(i, j) + temp * a(i, l);
					}
				}
			}
		}
		/*
		Form  C := alpha*A'*B  + beta*C
		*/
		else
		{
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < m; i++) {
					double temp = 0.0;
					for (int l = 0; l < k; l++) {
						temp += a[l][i] * b[l][j];
					}
					if (beta == 0.0) {
						c[i][j] = alpha * temp;
					}
					else {
						c[i][j] = alpha * temp + beta * c[i][j];
					}
				}
			}
		}
	}
	else
	{
		/*
		Form  C := alpha*A*B + beta*C.
		*/
		if (nota)
		{
			for (int j = 0; j < n; j++) {
				if (beta == 0) {
					for (int i = 0; i < m; i++) {
						c(i, j) = 0;
					}
				}
				else if (beta != 1) {
					for (int i = 0; i < m; i++) {
						c(i, j) = beta * c(i, j);
					}
				}
				for (int l = 0; l < k; l++) {
					temp = alpha * b(j, l);
					for (int i = 0; i < m; i++) {
						c(i, j) = c(i, j) + temp * a(i, l);
					}
				}
			}
		}

		/*
		Form  C := alpha*A'*B + beta*C
		*/
		else
		{
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < m; i++) {
					double temp = 0.0;
					for (int l = 0; l < k; l++) {
						temp += a[l][i] * b[j][l];
					}
					if (beta == 0.0) {
						c[i][j] = alpha * temp;
					}
					else {
						c[i][j] = alpha * temp + beta * c[i][j];
					}
				}
			}
		}
	}

	return;
}

//
// DGEMMTR  performs one of the matrix-matrix operations
//
//    C := alpha*op( A )*op( B ) + beta*C,
//
// where  op( X ) is one of
//
//    op( X ) = X   or   op( X ) = X**T,
//
// alpha and beta are scalars, and A, B and C are matrices, with op( A )
// an n by k matrix,  op( B )  a  k by n matrix and  C an n by n matrix.
// Thereby, the routine only accesses and updates the upper or lower
// triangular part of the result matrix C. This behaviour can be used if
// the resulting matrix C is known to be symmetric.
//
void cblas_dgemmtr(char uplo, char transa, char transb,
	int n, int k, double alpha, double* pA, int lda,
	double* pB, int ldb, double beta,
	double* pC, int ldc)
{
	double temp;
	int i, info, j, l, nrowa, nrowb, istart, istop;
	bool nota, notb, upper;
	double one = 1.0, zero = 0.0;
	//
	//     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
	//     transposed and set  NROWA and NROWB  as the number of rows of  A
	//     and  B  respectively.
	//
	nota = cblas_lsame(transa, 'N');
	notb = cblas_lsame(transb, 'N');
	if (nota) {
		nrowa = n;
	}
	else {
		nrowa = k;
	}
	if (notb) {
		nrowb = k;
	}
	else {
		nrowb = n;
	}
	upper = cblas_lsame(uplo, 'U');
	//
	//     Test the input parameters.
	//
	info = 0;
	if ((!upper) && (!cblas_lsame(uplo, 'L'))) {
		info = 1;
	}
	else if ((!nota) && (!cblas_lsame(transa, 'C')) &&
		(!cblas_lsame(transa, 'T'))) {
		info = 2;
	}
	else if ((!notb) && (!cblas_lsame(transb, 'C')) &&
		(!cblas_lsame(transb, 'T'))) {
		info = 3;
	}
	else if (n < 0) {
		info = 4;
	}
	else if (k < 0) {
		info = 5;
	}
	else if (lda < std::max(1, nrowa)) {
		info = 8;
	}
	else if (ldb < std::max(1, nrowb)) {
		info = 10;
	}
	else if (ldc < std::max(1, n)) {
		info = 13;
	}
	if (info != 0) {
		cblas_xerbla("DGEMMTR", info);
		return;
	}
	//
	//     Quick return if possible.
	//
	if (n == 0)
	{
		return;
	}

	MData a(nrowa, nrowa, pA, lda);
	MData b(nrowb, nrowb, pB, ldb);
	MData c(n, n, pC, ldc);

	//
	//     And if  alpha.eq.zero.
	//
	if (alpha == zero) {
		if (beta == zero) {
			for (j = 0; j < n; j++) {
				if (upper) {
					istart = 0;
					istop = j;
				}
				else {
					istart = j;
					istop = n;
				}
				for (i = istart; i < istop; i++) {
					c(i, j) = zero;
				}
			}
		}
		else {
			for (j = 0; j < n; j++) {
				if (upper) {
					istart = 0;
					istop = j;
				}
				else {
					istart = j;
					istop = n;
				}

				for (i = istart; i < istop; i++) {
					c(i, j) = beta * c(i, j);
				}
			}
		}
		return;
	}
	//
	//     Start the operations.
	//
	if (notb) {
		if (nota) {
			//
			//  Form  C := alpha*A*B + beta*C.
			//
			for (j = 0; j < n; j++) {
				if (upper) {
					istart = 0;
					istop = j + 1;
				}
				else {
					istart = j;
					istop = n;
				}
				if (beta == zero) {
					for (i = istart; i < istop; i++) {
						c(i, j) = zero;
					}
				}
				else if (beta != one) {
					for (i = istart; i < istop; i++) {
						c(i, j) = beta * c(i, j);
					}
				}
				for (l = 0; l < k; l++) {
					temp = alpha * b(l, j);
					for (i = istart; i < istop; i++) {
						c(i, j) = c(i, j) + temp * a(i, l);  
					}
				}
			}
		}
		else {
			//
			//   Form  C := alpha*A**T*B + beta*C
			//
			for (j = 0; j < n; j++) {
				if (upper) {
					istart = 0;
					istop = j + 1;
				}
				else {
					istart = j;
					istop = n;
				}

				for (i = istart; i < istop; i++) {
					temp = zero;
					for (l = 0; l < k; l++) {
						temp = temp + a(l, i) * b(l, j);
					}
					if (beta == zero) {
						c(i, j) = alpha * temp;
					}
					else {
						c(i, j) = alpha * temp + beta * c(i, j);
					}
				}
			}
		}
	}
	else {
		if (nota) {
			//
			//  Form  C := alpha*A*B**T + beta*C
			//
			for (j = 0; j < n; j++) {
				if (upper) {
					istart = 0;
					istop = j + 1;
				}
				else {
					istart = j;
					istop = n;
				}

				if (beta == zero) {
					for (i = istart; i < istop; i++) {
						c(i, j) = zero;
					}
				}
				else if (beta != one) {
					for (i = istart; i < istop; i++) {
						c(i, j) = beta * c(i, j);
					}
				}
				for (l = 0; l < k; l++) {
					temp = alpha * b(j, l);
					for (i = istart; i < istop; i++) {
						c(i, j) = c(i, j) + temp * a(i, l);
					}
				}
			}
		}
		else {
			//
			//           Form  C := alpha*A**T*B**T + beta*C
			//
			for (j = 0; j < n; j++) {
				if (upper) {
					istart = 0;
					istop = j + 1;
				}
				else {
					istart = j;
					istop = n;
				}

				for (i = istart; i < istop; i++) {
					temp = zero;
					for (l = 0; l < k; l++) {
						temp = temp + a(l, i) * b(j, l);
					}
					if (beta == zero) {
						c(i, j) = alpha * temp;
					}
					else {
						c(i, j) = alpha * temp + beta * c(i, j);
					}
				}
			}
		}
	}
	// 
	return;
}
/******************************************************************************/



/******************************************************************************/
/*
  Purpose:

    DTRMM performs B:=A*B or B:=B*A, A triangular, B rectangular.

  Discussion:

    This routine performs one of the matrix-matrix operations
      B := alpha*op( A )*B,
    or
      B := alpha*B*op( A ),
    where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
    non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
      op( A ) = A   or   op( A ) = A'.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 May 2024

  Author:

    This C version translated by Yang from LAPACK 3.12.0.

  Parameters:

    Input, char SIDE, specifies whether op(A) multiplies B from
    the left or right as follows:
    'L' or 'l': B := alpha*op( A )*B.
    'R' or 'r': B := alpha*B*op( A ).

    Input, char UPLO, specifies whether the matrix A is an upper or
    lower triangular matrix as follows:
    'U' or 'u': A is an upper triangular matrix.
    'L' or 'l': A is a lower triangular matrix.

    Input, char TRANS, specifies the form of op( A ) to be used in
    the matrix multiplication as follows:
    'N' or 'n': op( A ) = A.
    'T' or 't': op( A ) = A'.
    'C' or 'c': op( A ) = A'.

    Input, char DIAG, specifies whether or not A is unit triangular
    as follows:
    'U' or 'u': A is assumed to be unit triangular.
    'N' or 'n': A is not assumed to be unit triangular.

    Input, int M, the number of rows of B.  0 <= M.

    Input, int N, the number of columns of B.  
    0 <= N.

    Input, double ALPHA, the scalar  alpha.  When alpha is
    0.0, A is not referenced and B need not be set before entry.

    Input, double A[LDA*K], where k is m when  SIDE = 'L' or 'l'  
    and is  n  when  SIDE = 'R' or 'r'.
    Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
    upper triangular part of the array  A must contain the upper
    triangular matrix  and the strictly lower triangular part of
    A is not referenced.
    Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
    lower triangular part of the array  A must contain the lower
    triangular matrix  and the strictly upper triangular part of
    A is not referenced.
    Note that when  DIAG = 'U' or 'u',  the diagonal elements of
    A  are not referenced either,  but are assumed to be  unity.

    Input, integer LDA, the first dimension of A as declared
    in the calling program.  When SIDE = 'L' or 'l' then LDA must be at 
    least max ( 1, M ); when SIDE = 'R' or 'r', LDA must be at least 
    max ( 1, N ).

    Input/output, double B[LDB*N].
    Before entry, the leading m by n part of the array  B must contain 
    the matrix  B, and on exit is overwritten by the transformed matrix.

    Input, integer LDB, the first dimension of B as declared
    in the calling program.   max ( 1, M ) <= LDB.

	Test Cases:

A:

1.000 12.000 13.000 14.000
21.000 1.000 3.000 4.000
31.000 32.000 1.000 4.000
41.000 42.000 43.000 1.000

B:

1.000 1.000 1.000 1.000
2.000 1.000 3.000 4.000
3.000 2.000 1.000 4.000
4.000 4.000 3.000 1.000

	char side = 'L';
	char uplo = 'U';
	char transa = 'N';
	char diag = 'U';
	int m = 4;
	int n = 4;
	double alpha = 1.0;

After trmm:

120.000 95.000 92.000 115.000
27.000 23.000 18.000 20.000
19.000 18.000 13.000 8.000
4.000 4.000 3.000 1.000


	char side = 'L';
	char uplo = 'U';
	char transa = 'T';
	char diag = 'U';
	int m = 4;
	int n = 4;
	double alpha = 1.0;

After trmm:

1.000 1.000 1.000 1.000
14.000 13.000 15.000 16.000
22.000 18.000 23.000 29.000
38.000 30.000 33.000 47.000


	char side = 'R';
	char uplo = 'U';
	char transa = 'T';
	char diag = 'U';
	int m = 4;
	int n = 4;
	double alpha = 1.0;

After trmm:

40.000 8.000 5.000 1.000
109.000 26.000 19.000 4.000
96.000 21.000 17.000 4.000
105.000 17.000 7.000 1.000


	char side = 'R';
	char uplo = 'L';
	char transa = 'T';
	char diag = 'U';
	int m = 4;
	int n = 4;
	double alpha = 1.0;

After trmm:

1.000 22.000 64.000 127.000
2.000 43.000 97.000 257.000
3.000 65.000 158.000 254.000
4.000 88.000 255.000 462.000

*/
void cblas_dtrmm(char side, char uplo, char transa, char diag, int m, int n,
	double alpha, double* pA, int lda, double* pB, int ldb)
{
	int i;
	int info;
	int j;
	int k;
	int lside;
	int nounit;
	int nrowa;
	double temp = 0.0;
	int upper;
	/*
	  Test the input parameters.
	*/
	lside = cblas_lsame(side, 'L');

	if (lside)
	{
		nrowa = m;
	}
	else
	{
		nrowa = n;
	}

	nounit = cblas_lsame(diag, 'N');
	upper = cblas_lsame(uplo, 'U');

	info = 0;
	if (!lside && !cblas_lsame(side, 'R'))
	{
		info = 1;
	}
	else if (!upper && !cblas_lsame(uplo, 'L'))
	{
		info = 2;
	}
	else if (!cblas_lsame(transa, 'N') && !cblas_lsame(transa, 'T') &&
		!cblas_lsame(transa, 'C'))
	{
		info = 3;
	}
	else if (!cblas_lsame(diag, 'U') && !cblas_lsame(diag, 'N'))
	{
		info = 4;
	}
	else if (m < 0)
	{
		info = 5;
	}
	else if (n < 0)
	{
		info = 6;
	}
	else if (lda < std::max(1, nrowa))
	{
		info = 9;
	}
	else if (ldb < std::max(1, m))
	{
		info = 11;
	}

	if (info != 0)
	{
		cblas_xerbla("DTRMM", info);
		return;
	}
	/*
	  Quick return if possible.
	*/
	if (n == 0)
	{
		return;
	}

	MData a(nrowa, lda, pA, lda);
	//a.setData(pA);

	MData b(m, ldb, pB, ldb);
	//b.setData(pB);

	/*
	  And when alpha is 0.0.
	*/
	if (alpha == 0.0)
	{
		for (j = 0; j < n; j++)
		{
			for (i = 0; i < m; i++)
			{
				b[i][j] = 0.0;
			}
		}
		return;
	}
	/*
	  Start the operations.
	*/
	if (lside)
	{
		/*
		  Form  B := alpha*A*B.
		*/
		if (cblas_lsame(transa, 'N'))
		{
			if (upper)
			{
				/**

				|A A A|  B
				|  A A|  B
				|    A|  B
				*/

				for (j = 0; j < n; j++)
				{
					for (k = 0; k < m; k++)
					{
						if (b[k][j] != 0.0)
						{
							temp = alpha * b[k][j];
							for (i = 0; i < k; i++)
							{
								b[i][j] = b[i][j] + temp * a[i][k];
							}
							if (nounit)
							{
								temp = temp * a[k][k];
							}
							b[k][j] = temp;
						}
					}
				}
			}
			else
			{
				/**

				|A    |  B
				|A A  |  B
				|A A A|  B
				*/

				for (j = 0; j < n; j++)
				{
					for (k = m - 1; 0 <= k; k--)
					{
						if (b[k][j] != 0.0)
						{
							temp = alpha * b[k][j];
							b[k][j] = temp;
							if (nounit)
							{
								b[k][j] = b[k][j] * a[k][k];
							}
							for (i = k + 1; i < m; i++)
							{
								b[i][j] = b[i][j] + temp * a[i][k];
							}
						}
					}
				}

			}
		}
		else
		{
			/*
			  Form  B := alpha*A'*B.
			*/
			if (upper)
			{
				for (j = 0; j < n; j++)
				{
					for (i = m - 1; i >= 0; i--)
					{
						temp = b[i][j];
						if (nounit)
						{
							temp = temp * a[i][i];
						}
						for (k = 0; k < i; k++)
						{
							temp = temp + a[k][i] * b[k][j];
						}
						b[i][j] = alpha * temp;
					}
				}
			}
			else
			{
				for (j = 0; j < n; j++)
				{
					for (i = 0; i < m; i++)
					{
						temp = b[i][j];
						if (nounit)
						{
							temp = temp * a[i][i];
						}
						for (k = i + 1; k < m; k++)
						{
							temp = temp + a[k][i] * b[k][j];
						}
						b[i][j] = alpha * temp;
					}
				}
			}
		}
	}
	else
	{
		if (cblas_lsame(transa, 'n'))
		{
			/*
				Form  B := alpha*B*A.
			*/
			if (upper)
			{

				for (j = n - 1; j >= 0; j--)
				{
					temp = alpha;
					if (nounit)
					{
						temp = temp * a[j][j];
					}
					for (i = 0; i < m; i++)
					{
						b[i][j] = temp * b[i][j];
					}
					for (k = 0; k < j; k++)
					{
						if (a[k][j] != 0.0)
						{
							temp = alpha * a[k][j];
						}
						for (i = 0; i < m; i++)
						{
							b[i][j] = b[i][j] + temp * b[i][k];
						}
					}
				}

			}
			else
			{
				for (j = 0; j < n; j++)
				{
					temp = alpha;
					if (nounit)
					{
						temp = temp * a[j][j];
					}

					for (i = 0; i < m; i++)
					{
						b[i][j] = temp * b[i][j];
					}
					for (k = j + 1; k < n; k++)
					{
						if (a[k][j] != 0.0)
						{
							temp = alpha * a[k][j];
						}
						for (i = 0; i < m; i++)
						{
							b[i][j] = b[i][j] + temp * b[i][k];
						}
					}
				}

			}
		}
		else
		{
			/*
			  Form  B := alpha*B*A'.
			*/
			if (upper)
			{

				for (k = 0; k < n; k++)
				{
					for (j = 0; j < k; j++)
					{
						if (a[j][k] != 0.0)
						{
							temp = alpha * a[j][k];
							for (i = 0; i < m; i++)
							{
								b[i][j] = b[i][j] + temp * b[i][k];
							}
						}
					}
					temp = alpha;
					if (nounit)
					{
						temp = temp * a[k][k];
					}
					if (temp != 1.0)
					{
						for (i = 0; i < m; i++)
						{
							b[i][k] = temp * b[i][k];
						}
					}
				}

			}
			else
			{
				for (k = n - 1; k >= 0; k--)
				{
					for (j = k + 1; j < n; j++)
					{
						if (a[j][k] != 0.0)
						{
							temp = alpha * a[j][k];
							for (i = 0; i < m; i++)
							{
								b[i][j] = b[i][j] + temp * b[i][k];
							}
						}
					}
					temp = alpha;
					if (nounit)
					{
						temp = temp * a[k][k];
					}
					if (temp != 1.0)
					{
						for (i = 0; i < m; i++)
						{
							b[i][k] = temp * b[i][k];
						}
					}
				}

			}
		}
	}
	return;
}
/******************************************************************************/



/******************************************************************************/
/*
  Purpose:

	DTRSM solves A*X=alpha*B or X*A=alpha*B, for triangular A, rectangular B.

  Discussion:

	DTRSM solves one of the matrix equations
	  op( A )*X = alpha*B,
	or
	  X*op( A ) = alpha*B,
	where alpha is a scalar, X and B are m by n matrices, A is a unit, or
	non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
	  op( A ) = A   or   op( A ) = A'.
	The matrix X is overwritten on B.

  Licensing:

	This code is distributed under the GNU LGPL license.

  Modified:

	06 April 2014

  Author:

	This C version by John Burkardt.

  Parameters:

	Input, char SIDE, specifies whether op( A ) appears on the left
	or right of X as follows:
	'L' or 'l': op( A )*X = alpha*B. -- the row of op(A) shoule be m

	'R' or 'r': X*op( A ) = alpha*B. -- the row of op(A) shoule be n

	Input, char UPLO, specifies whether the matrix A is an upper or
	lower triangular matrix as follows:
	'U' or 'u': A is an upper triangular matrix.
	'L' or 'l': A is a lower triangular matrix.

	Input, char TRANSA, specifies the form of op( A ) to be used in
	the matrix multiplication as follows:
	'N' or 'n': op( A ) = A.
	'T' or 't': op( A ) = A'.
	'C' or 'c': op( A ) = A'.

	Input, char DIAG, specifies whether or not A is unit triangular
	as follows:
	'U' or 'u': A is assumed to be unit triangular.
	'N' or 'n': A is not assumed to be unit triangular.

	Input, int M, the number of rows of B.  0 <= M.

	Input, int N, the number of columns of B.  0 <= N.

	Input, double ALPHA, the scalar alpha.  When alpha is
	0.0 then A is not referenced and B need not be set before entry.

	Input, double A[LDA*K] where K is M when SIDE = 'L' or 'l'
	and K is N when SIDE = 'R' or 'r'.
	Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
	upper triangular part of the array  A must contain the upper
	triangular matrix  and the strictly lower triangular part of
	A is not referenced.
	Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
	lower triangular part of the array  A must contain the lower
	triangular matrix  and the strictly upper triangular part of
	A is not referenced.
	Note that when  DIAG = 'U' or 'u',  the diagonal elements of
	A  are not referenced either,  but are assumed to be  unity.

	Input, int LDA, the first dimension of A as declared
	in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
	LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
	then LDA must be at least max( 1, n ).

	Input/output, double B[M*LDB].
	Before entry, B is a m by n matrix. the leading m by n part of the array B must
	contain the right-hand side matrix B, and on exit is
	overwritten by the solution matrix X.

	Input, int LDB, the first dimension of B as declared
	in the calling program.  LDB must be at least max ( 1, N ).
*/
void cblas_dtrsm(char side, char uplo, char transa, char diag, int m, int n,
	double alpha, double* pA, int lda, double* pB, int ldb)
{
	int i;
	int info;
	int j;
	int k;
	int lside;
	int nounit;
	int nrowa;
	int ncola;
	double temp;
	int upper;
	/*
	  Test the input parameters.
	*/
	lside = cblas_lsame(side, 'L');

	if (lside)
	{
		nrowa = m;
		ncola = n;
	}
	else
	{
		nrowa = n;
		ncola = m;
	}

	//
	printf("  nrowa: %d  \n", nrowa);
	printf("  ncola: %d  \n", ncola);
	//

	// whether or not A is unit triangular
	nounit = cblas_lsame(diag, 'N');

	// whether the matrix A is an upper or
	// lower triangular matrix
	upper = cblas_lsame(uplo, 'U');

	info = 0;

	if ((!lside) && (!cblas_lsame(side, 'R')))
	{
		info = 1;
	}
	else if ((!upper) && (!cblas_lsame(uplo, 'L')))
	{
		info = 2;
	}
	else if ((!cblas_lsame(transa, 'N')) &&
		(!cblas_lsame(transa, 'T')) &&
		(!cblas_lsame(transa, 'C')))
	{
		info = 3;
	}
	else if ((!cblas_lsame(diag, 'U')) && (!cblas_lsame(diag, 'N')))
	{
		info = 4;
	}
	else if (m < 0)
	{
		info = 5;
	}
	else if (n < 0)
	{
		info = 6;
	}
	// C programming language uses column-major storing style 
	else if (lda < std::max(1, ncola))
	{
		info = 9;
	}
	else if (ldb < std::max(1, n))
	{
		info = 11;
	}

	if (info != 0)
	{
		cblas_xerbla("DTRSM", info);
		return;
	}
	/*
	  Quick return if possible.
	*/
	if (n == 0)
	{
		return;
	}

	MData A;
	if (upper)
	{
		A = MData(nrowa, lda, pA, lda, 4);
	}
	else {
		A = MData(nrowa, lda, pA, lda, 5);
	}
	A.show();

	MData a(nrowa, lda, pA, lda);
	MData b(m, ldb, pB, ldb);

	//printf("    A:    \n" );
	//a.show();

	/*
	  and when alpha is 0.0.
	*/
	if (alpha == 0.0)
	{
		for (j = 0; j < n; j++)
		{
			for (i = 0; i < m; i++)
			{
				b[i][j] = 0.0;
			}
		}
		return;
	}
	/*
	  Start the operations.
	*/
	if (lside)
	{
		/*
		  Form  B := alpha*inv( a )*B.
		*/
		if (cblas_lsame(transa, 'N'))
		{
			//  Debug -- BEGIN -- 
			printf("     Trans: %c  \n", transa);
			printf("     Form  B := alpha*inv( a )*B. \n");
			//  Debug -- END --  

			if (upper)
			{
				for (j = 0; j < n; j++)
				{
					if (alpha != 1.0)
					{
						for (i = 0; i < m; i++)
						{
							b[i][j] = alpha * b[i][j];
						}
					}
					for (k = m - 1; 0 <= k; k--)
					{
						if (b[k][j] != 0.0)
						{
							if (nounit)
							{
								b[k][j] = b[k][j] / A(k, k);
							}
							for (i = 0; i < k; i++)
							{
								b[i][j] = b[i][j] - b[k][j] * A(i, k);
							}
						}
					}
				}
			}
			else
			{
				for (j = 0; j < n; j++)
				{
					if (alpha != 1.0)
					{
						for (i = 0; i < m; i++)
						{
							b[i][j] = alpha * b[i][j];
						}
					}
					for (k = 0; k < m; k++)
					{
						if (b[k][j] != 0.0)
						{
							if (nounit)
							{
								b[k][j] = b[k][j] / A(k, k);
							}
							for (i = k + 1; i < m; i++)
							{
								b[i][j] = b[i][j] - b[k][j] * A(i, k);
							}
						}
					}
				}
			}
		}
		/*
		  Form  B := alpha*inv( A' )*B.
		*/
		else
		{
			//  Debug -- BEGIN -- 
			printf("     Trans: %c  \n", transa);
			printf("     Form  B := alpha*inv( A' )*B. \n");
			//  Debug -- END --  
			if (upper)
			{
				for (j = 0; j < n; j++)
				{
					for (i = 0; i < m; i++)
					{
						temp = alpha * b[i][j];
						for (k = 0; k < i; k++)
						{
							temp = temp - A(k, i) * b[k][j];
						}
						if (nounit)
						{
							temp = temp / A(i, i);
						}
						b[i][j] = temp;
					}
				}
			}
			else
			{
				for (j = 0; j < n; j++)
				{
					for (i = m - 1; 0 <= i; i--)
					{
						temp = alpha * b[i][j];
						for (k = i + 1; k < m; k++)
						{
							temp = temp - A(k, i) * b[k][j];
						}
						if (nounit)
						{
							temp = temp / A(i, i);
						}
						b[i][j] = temp;
					}
				}
			}
		}
	}
	/*
	  Form  B := alpha*B*inv( A ).
	*/
	else
	{
		if (cblas_lsame(transa, 'N'))
		{
			//  Debug -- BEGIN -- 
			printf("     Trans: %c  \n", transa);
			printf("     Form  B := alpha*B*inv( A ). \n");
			//  Debug -- END --  

			if (upper)
			{
				//  Debug -- BEGIN -- 
				printf("     Upper: U -- %c  \n", transa);
				//  Debug -- END --  

				for (j = 0; j < n; j++)
				{
					if (alpha != 1.0)
					{
						for (i = 0; i < m; i++)
						{
							b[i][j] = alpha * b[i][j];
						}
					}
					for (k = 0; k < j; k++)
					{
						if (A(k, j) != 0.0)
						{
							for (i = 0; i < m; i++)
							{
								b[i][j] = b[i][j] - A(k, j) * b[i][k];
							}
						}
					}
					if (nounit)
					{
						temp = 1.0 / A(j, j);
						for (i = 0; i < m; i++)
						{
							b[i][j] = temp * b[i][j];
						}
					}
				}
			}
			else
			{
				//  Debug -- BEGIN -- 
				printf("     Upper: L -- %c  \n", transa);
				//  Debug -- END --  
				for (j = n - 1; 0 <= j; j--)
				{
					if (alpha != 1.0)
					{
						for (i = 0; i < m; i++)
						{
							b[i][j] = alpha * b[i][j];
						}
					}
					for (k = j + 1; k < n; k++)
					{
						if (A(k, j) != 0.0)
						{
							for (i = 0; i < m; i++)
							{
								b[i][j] = b[i][j] - A(k, j) * b[i][k];
							}
						}
					}
					if (nounit)
					{
						temp = 1.0 / A(j, j);
						for (i = 0; i < m; i++)
						{
							b[i][j] = temp * b[i][j];
						}
					}
				}
			}
		}
		/*
		  Form  B := alpha*B*inv( A' ).
		*/
		else
		{
			//  Debug -- BEGIN -- 
			printf("     Trans: %c  \n", transa);
			printf("     Form  B := alpha*B*inv( A' ). \n");
			//  Debug -- END --  

			if (upper)
			{
				for (k = n - 1; 0 <= k; k--)
				{
					if (nounit)
					{
						temp = 1.0 / A(k, k);
						for (i = 0; i < m; i++)
						{
							b[i][k] = temp * b[i][k];
						}
					}
					for (j = 0; j < k; j++)
					{
						if (A(j, k) != 0.0)
						{
							temp = A(j, k);
							for (i = 0; i < m; i++)
							{
								b[i][j] = b[i][j] - temp * b[i][k];
							}
						}
					}
					if (alpha != 1.0)
					{
						for (i = 0; i < m; i++)
						{
							b[i][k] = alpha * b[i][k];
						}
					}
				}
			}
			else
			{
				for (k = 0; k < n; k++)
				{
					if (nounit)
					{
						temp = 1.0 / A(k, k);
						for (i = 0; i < m; i++)
						{
							b[i][k] = temp * b[i][k];
						}
					}
					for (j = k + 1; j < n; j++)
					{
						if (A(j, k) != 0.0)
						{
							temp = A(j, k);
							for (i = 0; i < m; i++)
							{
								b[i][j] = b[i][j] - temp * b[i][k];
							}
						}
					}
					if (alpha != 1.0)
					{
						for (i = 0; i < m; i++)
						{
							b[i][k] = alpha * b[i][k];
						}
					}
				}
			}
		}
	}
	return;
}





/**
@brief
DSYRK  performs one of the symmetric rank k operations
C := alpha*A*A^T + beta*C,
or
C := alpha*A^T*A + beta*C,
where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
and  A  is an  n by k  matrix in the first case and a  k by n  matrix
in the second case.
https://netlib.org/lapack/explore-html/dc/d05/dsyrk_8f_source.html

@param[in] UPLO char
On  entry,   UPLO  specifies  whether  the  upper  or  lower
triangular  part  of the  array  C  is to be  referenced  as
follows:
UPLO = 'U' or 'u'   Only the  upper triangular part of  C
is to be referenced.

UPLO = 'L' or 'l'   Only the  lower triangular part of  C
is to be referenced.

@param[in] TRANS: char
On entry,  TRANS  specifies the operation to be performed as follows:
TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C.

@param[in]  N: int, specifies the order of the matrix C.  N must be
at least zero.
*/
void cblas_dsyrk(char uplo, char trans, int n, int k,
	double alpha, double *A, int lda, double beta, double *C, int ldc)
{
	double temp;
	int i, info, j, l, nrowa;
	bool upper;

	double ONE = 1.0;
	double ZERO = 0.0;

	/*
	*     Test the input parameters.
	*/
	if (cblas_lsame(trans, 'N')) {
		nrowa = n;
	}
	else {
		nrowa = k;
	}
	upper = cblas_lsame(uplo, 'U');

	info = 0;
	if ((!upper) & (!cblas_lsame(uplo, 'L'))) {
		info = 0;
	}
	else if ((!cblas_lsame(trans, 'N'))  &
		(!cblas_lsame(trans, 'T'))  &
		(!cblas_lsame(trans, 'C'))) {
		info = 2;
	}
	else  if (n < 0) {
		info = 3;
	}
	else  if (k < 0) {
		info = 4;
	}
	else  if (lda < std::max(1, nrowa)) {
		info = 7;
	}
	else  if (ldc < std::max(1, n)) {
		info = 10;
	}
	if (info != 0) {
		xerbla("DSYRK ", info);
		return;
	}

	/*
	*     Quick return if possible.
	*/
	if ((n == 0) || (((alpha == ZERO) || (k == 0)) &  (beta == ONE))) {
		return;
	}

	/*
	*     And when alpha.eq.ZERO.
	*/
	if (alpha == ZERO) {
		if (upper) {
			if (beta == ZERO) {
				for (j = 0; j < n; j++) {
					for (i = 0; i < j; i++) {
						C[i * ldc + j] = ZERO;
					}
				}
			}
			else {
				for (j = 0; j < n; j++) {
					for (i = 0; i < j; i++) {
						C[i * ldc + j] = beta*C[i * ldc + j];
					}
				}
			}
		}
		else {
			if (beta == ZERO) {
				for (j = 0; j < n; j++) {
					for (i = j; i < n; i++) {
						C[i * ldc + j] = ZERO;
					}
				}
			}
			else {
				for (j = 0; j < n; j++) {
					for (i = j; i < n; i++) {
						C[i * ldc + j] = beta*C[i * ldc + j];
					}
				}
			}
		}
		return;
	}

	/*
	*     Start the operations.
	*/
	if (cblas_lsame(trans, 'N')) {
		// printf("  Form  C := alpha*A*A**T + beta*C. \n");

		/*
		*  Form  C := alpha*A*A**T + beta*C.
		*/
		if (upper) {

			// printf("    UPLO := U  \n");

			for (j = 0; j < n; j++) {

				if (beta == ZERO) {
					for (i = 0; i < j; i++) {
						C[i * ldc + j] = ZERO;
					}
				}
				else  if (beta != ONE) {
					for (i = 0; i < j; i++) {
						C[i * ldc + j] = beta*C[i * ldc + j];
					}
				}

				for (l = 0; l < k; l++) {
					if (A[j * lda + l] != ZERO) {
						temp = alpha*A[j * lda + l];
						for (i = 0; i <= j; i++) {
							C[i * ldc + j] = C[i * ldc + j] + temp*A[i * lda + l];
						}
					}
				}
			}
		}
		else {

			//printf("    UPLO := L  \n");
			 
			for (j = 0; j < n; j++) {
				if (beta == ZERO) {
					for (i = j; i < n; i++) {
						C[i * ldc + j] = ZERO;
					}
				}
				else  if (beta != ONE) {
					for (i = j; i < n; i++) {
						C[i * ldc + j] = beta*C[i * ldc + j];
					}
				}
				for (l = 0; l < k; l++) {
					if (A[j * lda + l] != ZERO) {
						temp = alpha*A[j * lda + l];
						for (i = j; i < n; i++) {
							C[i * ldc + j] = C[i * ldc + j] + temp*A[i * lda + l];
						}
					}
				}
			}
		}
	}
	else {
		//printf("  Form  C := alpha*A**T*A + beta*C. \n"); 
		/*
		*  Form  C := alpha*A**T*A + beta*C.
		*/
		if (upper) {
			for (j = 0; j < n; j++) {
				for (i = 0; i <= j; i++) {
					temp = ZERO;
					for (l = 0; l < k; l++) {
						temp = temp + A[l * lda + i] * A[l * lda + j];
					}
					if (beta == ZERO) {
						C[i * ldc + j] = alpha*temp;
					}
					else {
						C[i * ldc + j] = alpha*temp + beta*C[i * ldc + j];
					}
				}
			}
		}
		else {
			for (j = 0; j < n; j++) {
				for (i = j; i < n; i++) {
					temp = ZERO;
					for (l = 0; l < k; l++) {
						temp = temp + A[l * lda + i] * A[l * lda + j];
					}
					if (beta == ZERO) {
						C[i * ldc + j] = alpha*temp;
					}
					else {
						C[i * ldc + j] = alpha*temp + beta*C[i * ldc + j];
					}
				}
			}
		}
	}
}


