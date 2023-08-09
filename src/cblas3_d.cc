# include <stdio.h>
# include <stdlib.h>
# include <time.h>
#include "cblas.h"

#define min(x,y) (((x) < (y)) ? (x) : (y))
#define max(x,y) (((x) > (y)) ? (x) : (y)) 

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
	double alpha, double a[], int lda, double b[], int ldb, double beta,
	double c[], int ldc)
{
  int i;
  int j;
  int l;
  int nrowa;
  int nrowb;
  int nota;
  int notb;
  double temp;
/*
  Set NOTA and NOTB as true if A and B respectively are not
  transposed and set NROWA, NCOLA and NROWB as the number of rows
  and columns of A and the number of rows of B respectively.
*/
  nota = ( ( transa == 'N' ) || ( transa == 'n' ) );

  if ( nota )
  {
    nrowa = m;
  }
  else
  {
    nrowa = k;
  }

  notb = ( ( transb == 'N' ) || ( transb == 'n' ) );

  if ( notb )
  {
    nrowb = k;
  }
  else
  {
    nrowb = n;
  }
/*
  Test the input parameters.
*/
  if ( ! ( transa == 'N' || transa == 'n' ||
           transa == 'C' || transa == 'c' ||
           transa == 'T' || transa == 't' ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DGEMM - Fatal error!\n" );
    fprintf ( stderr, "  Input TRANSA had illegal value.\n" );
    //exit ( 1 );
	return;
  }

  if ( ! ( transb == 'N' || transb == 'n' ||
           transb == 'C' || transb == 'c' ||
           transb == 'T' || transb == 't' ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DGEMM - Fatal error!\n" );
    fprintf ( stderr, "  Input TRANSB had illegal value.\n" );
	//exit ( 1 );
	return;
  }

  if ( m < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DGEMM - Fatal error!\n" );
    fprintf ( stderr, "  Input M had illegal value.\n" );
	//exit ( 1 );
	return;
  }

  if ( n < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DGEMM - Fatal error!\n" );
    fprintf ( stderr, "  Input N had illegal value.\n" );
	//exit ( 1 );
	return;
  }

  if ( k  < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DGEMM - Fatal error!\n" );
    fprintf ( stderr, "  Input K had illegal value.\n" );
	//exit ( 1 );
	return;
  }

  //if ( lda < i4_max ( 1, nrowa ) )
  if (lda < max(1, nrowa))
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DGEMM - Fatal error!\n" );
    fprintf ( stderr, "  Input LDA had illegal value.\n" );
	//exit ( 1 );
	return;
  }

  if ( ldb < max ( 1, nrowb ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DGEMM - Fatal error!\n" );
    fprintf ( stderr, "  Input LDB had illegal value.\n" );
	//exit ( 1 );
	return;
  }

  if (ldc < max(1, m))//if ( ldc < i4_max ( 1, m ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DGEMM - Fatal error!\n" );
    fprintf ( stderr, "  Input LDC had illegal value.\n" );
	//exit ( 1 );
	return;
  }
/*
  Quick return if possible.
*/
  if ( m == 0 )
  {
    return;
  }

  if ( n == 0 )
  {
    return;
  }

  if ( ( alpha == 0.0 || k == 0 ) && ( beta == 1.0 ) )
  {
    return;
  }
/*
  And if alpha is 0.0.
*/
  if ( alpha == 0.0 )
  {
    if ( beta == 0.0 )
    {
      for ( j = 0; j < n; j++ )
      {
        for ( i = 0; i < m; i++ )
        {
          c[i+j*ldc] = 0.0;
        }
      }
    }
    else
    {
      for ( j = 0; j < n; j++ )
      {
        for ( i = 0; i < m; i++ )
        {
          c[i+j*ldc] = beta * c[i+j*ldc];
        }
      }
    }
    return;
  }
/*
  Start the operations.
*/
  if ( notb )
  {
	  if (nota)
	  {
		  for (j = 0; j < n; j++)
		  {
			  for (i = 0; i < m; i++)
			  {
				  temp = 0.0;
				  for (l = 0; l < k; l++)
				  {
					  temp = temp + a[l + i*lda] * b[j + l*ldb];
				  }
				  if (beta == 0.0)
				  {
					  c[i*ldc + j] = alpha * temp;
				  }
				  else
				  {
					  c[i*ldc + j] = alpha * temp + beta * c[i*ldc + j];
				  }
			  }
		  }


	  }
	  /*
	  Form  C := alpha*A'*B' + beta*C
	  */
	  else
	  {
		  for (j = 0; j < n; j++)
		  {
			  if (beta == 0.0)
			  {
				  for (i = 0; i < m; i++)
				  {
					  c[i*ldc + j] = 0.0;
				  }
			  }
			  else if (beta != 1.0)
			  {
				  for (i = 0; i < m; i++)
				  {
					  c[i*ldc + j] = beta * c[i*ldc + j];
				  }
			  }

			  for (l = 0; l < k; l++)
			  {
				  if (b[j + l*ldb] != 0.0)
				  {
					  temp = alpha * b[j + l*ldb];
					  for (i = 0; i < m; i++)
					  {
						  c[i*ldc + j] = c[i*ldc + j] + temp * a[i + l*lda];
					  }
				  }
			  }
		  }
	  }
  }
/*
  Form  C := alpha*A*B' + beta*C
*/
  else
  {
	  /*
	  Form  C := alpha*A*B + beta*C.
	  */
	  if (nota)
	  {
		  for (j = 0; j < n; j++)
		  {
			  for (i = 0; i < m; i++)
			  {
				  temp = 0.0;
				  for (l = 0; l < k; l++)
				  {
					  temp = temp + a[l + i*lda] * b[l + j*ldb];
				  }

				  if (beta == 0.0)
				  {
					  c[i*ldc + j] = alpha * temp;
				  }
				  else
				  {
					  c[i*ldc + j] = alpha * temp + beta * c[i*ldc + j];
				  }
			  }
		  }
	  }
	  /*
	  Form  C := alpha*A'*B + beta*C
	  */
	  else
	  {
		  for (j = 0; j < n; j++)
		  {
			  if (beta == 0.0)
			  {
				  for (i = 0; i < m; i++)
				  {
					  c[i*ldc + j] = 0.0;
				  }
			  }
			  else if (beta != 1.0)
			  {
				  for (i = 0; i < m; i++)
				  {
					  c[i*ldc + j] = beta * c[i*ldc + j];
				  }
			  }

			  for (l = 0; l < k; l++)
			  {
				  if (b[l + j*ldb] != 0.0)
				  {
					  temp = alpha * b[l + j*ldb];
					  for (i = 0; i < m; i++)
					  {
						  c[i*ldc + j] = c[i*ldc + j] + temp * a[i + l*lda];
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

    05 April 2014

  Author:

    This C version by John Burkardt.

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
*/
void cblas_dtrmm(char side, char uplo, char transa, char diag, int m, int n,
	double alpha, double a[], int lda, double b[], int ldb)
{
  int i;
  int info;
  int j;
  int k;
  int lside;
  int nounit;
  int nrowa;
  double temp;
  int upper;
/*
  Test the input parameters.
*/
  lside = cblas_lsame ( side, 'L' );

  if ( lside )
  {
    nrowa = m;
  }
  else
  {
    nrowa = n;
  }

  nounit = cblas_lsame ( diag, 'N' );
  upper = cblas_lsame ( uplo, 'U' );

  info = 0;
  if ( ! lside && ! cblas_lsame ( side, 'R' ) )
  {
    info = 1;
  }
  else if ( ! upper && ! cblas_lsame ( uplo, 'L' ) )
  {
    info = 2;
  }
  else if ( ! cblas_lsame ( transa, 'N' ) && ! cblas_lsame ( transa, 'T' ) && 
            ! cblas_lsame ( transa, 'C' ) )
  {
    info = 3;
  }
  else if ( ! cblas_lsame ( diag, 'U' ) && ! cblas_lsame ( diag, 'N' ) )
  {
    info = 4;
  }
  else if ( m < 0 )
  {
    info = 5;
  }
  else if ( n < 0 )
  {
    info = 6;
  }
  else if ( lda < max ( 1, nrowa ) )
  {
    info = 9;
  }
  else if ( ldb < max( 1, m ) )
  {
    info = 11;
  }

  if ( info != 0 ) 
  {
    xerbla ( "DTRMM", info );
    return;
  }
/*
  Quick return if possible.
*/
  if ( n == 0 )
  {
    return;
  }
/*
  And when alpha is 0.0.
*/
  if ( alpha == 0.0 )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        b[i+j*ldb] = 0.0;
      }
    }
    return;
  }
/*
  Start the operations.
*/
  if ( lside )
  {
/*
  Form  B := alpha*A*B.
*/
    if ( cblas_lsame ( transa, 'N' ) )
    {
		if (upper)
		{
			for (j = 0; j < n; j++)
			{
				for (i = m - 1; 0 <= i; i--)
				{
					temp = b[i + j*ldb];
					if (nounit)
					{
						temp = temp * a[i + i*lda];
					}
					for (k = 0; k < i; k++)
					{
						temp = temp + a[k + i*lda] * b[k + j*ldb];
					}
					b[i + j*ldb] = alpha * temp;
				}
			}
		}
		else
		{
			for (j = 0; j < n; j++)
			{
				for (i = 0; i < m; i++)
				{
					temp = b[i + j*ldb];
					if (nounit)
					{
						temp = temp * a[i + i*lda];
					}
					for (k = i + 1; k < m; k++)
					{
						temp = temp + a[k + i*lda] * b[k + j*ldb];
					}
					b[i + j*ldb] = alpha * temp;
				}
			}
		}
    }
/*
  Form  B := alpha*A'*B.
*/
    else
    {
		if (upper)
		{
			for (j = 0; j < n; j++)
			{
				for (k = 0; k < m; k++)
				{
					if (b[k + j*ldb] != 0.0)
					{
						temp = alpha * b[k + j*ldb];
						for (i = 0; i < k; i++)
						{
							b[i + j*ldb] = b[i + j*ldb] + temp * a[i + k*lda];
						}
						if (nounit)
						{
							temp = temp * a[k + k*lda];
						}
						b[k + j*ldb] = temp;
					}
				}
			}
		}
		else
		{
			for (j = 0; j < n; j++)
			{
				for (k = m - 1; 0 <= k; k--)
				{
					if (b[k + j*ldb] != 0.0)
					{
						temp = alpha * b[k + j*ldb];
						b[k + j*ldb] = temp;
						if (nounit)
						{
							b[k + j*ldb] = b[k + j*ldb] * a[k + k*lda];
						}
						for (i = k + 1; i < m; i++)
						{
							b[i + j*ldb] = b[i + j*ldb] + temp * a[i + k*lda];
						}
					}
				}
			}
		}

    }
  }
/*
  Form  B := alpha*B*A.
*/
  else
  {
    if ( cblas_lsame ( transa, 'n' ) )
    {
		if (upper)
		{
			for (k = 0; k < n; k++)
			{
				for (j = 0; j < k; j++)
				{
					if (a[j + k*lda] != 0.0)
					{
						temp = alpha * a[j + k*lda];
						for (i = 0; i < m; i++)
						{
							b[i + j*ldb] = b[i + j*ldb] + temp * b[i + k*ldb];
						}
					}
				}
				temp = alpha;
				if (nounit)
				{
					temp = temp * a[k + k*lda];
				}
				if (temp != 1.0)
				{
					for (i = 0; i < m; i++)
					{
						b[i + k*ldb] = temp * b[i + k*ldb];
					}
				}
			}
		}
		else
		{
			for (k = n - 1; 0 <= k; k--)
			{
				for (j = k + 1; j < n; j++)
				{
					if (a[j + k*lda] != 0.0)
					{
						temp = alpha * a[j + k*lda];
						for (i = 0; i < m; i++)
						{
							b[i + j*ldb] = b[i + j*ldb] + temp * b[i + k*ldb];
						}
					}
				}
				temp = alpha;
				if (nounit)
				{
					temp = temp * a[k + k*lda];
				}
				if (temp != 1.0)
				{
					for (i = 0; i < m; i++)
					{
						b[i + k*ldb] = temp * b[i + k*ldb];
					}
				}
			}
		}
    }
/*
  Form  B := alpha*B*A'.
*/
    else
    {
		if (upper)
		{
			for (j = n - 1; 0 <= j; j--)
			{
				temp = alpha;
				if (nounit)
				{
					temp = temp * a[j + j*lda];
				}
				for (i = 0; i < m; i++)
				{
					b[i + j*ldb] = temp * b[i + j*ldb];
				}
				for (k = 0; k < j; k++)
				{
					if (a[k + j*lda] != 0.0)
					{
						temp = alpha * a[k + j*lda];
						for (i = 0; i < m; i++)
						{
							b[i + j*ldb] = b[i + j*ldb] + temp * b[i + k*ldb];
						}
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
					temp = temp * a[j + j*lda];
				}
				for (i = 0; i < m; i++)
				{
					b[i + j*ldb] = temp * b[i + j*ldb];
				}
				for (k = j + 1; k < n; k++)
				{
					if (a[k + j*lda] != 0.0)
					{
						temp = alpha * a[k + j*lda];
						for (i = 0; i < m; i++)
						{
							b[i + j*ldb] = b[i + j*ldb] + temp * b[i + k*ldb];
						}
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
    'L' or 'l': op( A )*X = alpha*B.
    'R' or 'r': X*op( A ) = alpha*B.

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

    Input/output, double B[LDB*N].
    Before entry, the leading m by n part of the array B must
    contain the right-hand side matrix B, and on exit is
    overwritten by the solution matrix X.

    Input, int LDB, the first dimension of B as declared
    in the calling program.  LDB must be at least max ( 1, M ).
*/
void cblas_dtrsm(char side, char uplo, char transa, char diag, int m, int n,
	double alpha, double a[], int lda, double b[], int ldb)
{
  int i;
  int info;
  int j;
  int k;
  int lside;
  int nounit;
  int nrowa;
  double temp;
  int upper;
/*
  Test the input parameters.
*/
  lside = cblas_lsame ( side, 'L' );

  if ( lside )
  {
    nrowa = m;
  }
  else
  {
    nrowa = n;
  }

  nounit = cblas_lsame ( diag, 'N' );
  upper = cblas_lsame ( uplo, 'U' );

  info = 0;

  if ( ( ! lside ) && ( ! cblas_lsame ( side, 'R' ) ) )
  {
    info = 1;
  }
  else if ( ( ! upper ) && ( ! cblas_lsame ( uplo, 'L' ) ) )
  {
    info = 2;
  }
  else if ( ( ! cblas_lsame ( transa, 'N' ) ) && 
            ( ! cblas_lsame ( transa, 'T' ) ) && 
            ( ! cblas_lsame ( transa, 'C' ) ) )
  {
    info = 3;
  }
  else if ( ( ! cblas_lsame ( diag, 'U' ) ) && ( ! cblas_lsame ( diag, 'N' ) ) )
  {
    info = 4;
  }
  else if ( m < 0 )
  {
    info = 5;
  }
  else if ( n < 0 )
  {
    info = 6;
  }
  else if ( lda < max( 1, nrowa ) )
  {
    info = 9;
  }
  else if ( ldb < max( 1, m ) )
  {
    info = 11;
  }

  if ( info != 0 ) 
  {
    xerbla ( "DTRSM", info );
    return;
  }
/*
  Quick return if possible.
*/
  if ( n == 0 )
  {
    return;
  }
/*
  and when alpha is 0.0.
*/
  if ( alpha == 0.0 )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        b[i+j*ldb] = 0.0;
      }
    }
    return;
  }
/*
  Start the operations.
*/
  if ( lside ) 
  {
/*
  Form  B := alpha*inv( a )*B.
*/
    if ( cblas_lsame ( transa, 'N' ) )
    {
      if ( upper )
      {
        for ( j = 0; j < n; j++ )
        {
          if ( alpha != 1.0 )
          {
            for ( i = 0; i < m; i++ )
            {
              b[i+j*ldb] = alpha * b[i+j*ldb];
            }
          }
          for ( k = m - 1; 0 <= k; k-- )
          {
            if ( b[k+j*ldb] != 0.0 )
            {
              if ( nounit )
              {
                b[k+j*ldb] = b[k+j*ldb] / a[k+k*lda];
              }
              for ( i = 0; i < k; i++ )
              {
                b[i+j*ldb] = b[i+j*ldb] - b[k+j*ldb] * a[i+k*lda];
              }
            }
          }
        }
      }
      else
      {
        for ( j = 0; j < n; j++ )
        {
          if ( alpha != 1.0 )
          {
            for ( i = 0; i < m; i++ )
            {
              b[i+j*ldb] = alpha * b[i+j*ldb];
            }
          }
          for ( k = 0; k < m; k++ )
          {
            if ( b[k+j*ldb] != 0.0 )
            {
              if ( nounit )
              {
                b[k+j*ldb] = b[k+j*ldb] / a[k+k*lda];
              }
              for ( i = k + 1; i < m; i++ )
              {
                b[i+j*ldb] = b[i+j*ldb] - b[k+j*ldb] * a[i+k*lda];
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
      if ( upper )
      {
        for ( j = 0; j < n; j++ )
        {
          for ( i = 0; i < m; i++ )
          {
            temp = alpha * b[i+j*ldb];
            for ( k = 0; k < i; k++ )
            {
              temp = temp - a[k+i*lda] * b[k+j*ldb];
            }
            if ( nounit )
            {
              temp = temp / a[i+i*lda];
            }
            b[i+j*ldb] = temp;
          }
        }
      }
      else
      {
        for ( j = 0; j < n; j++ )
        {
          for ( i = m - 1; 0 <= i; i-- )
          {
            temp = alpha * b[i+j*ldb];
            for ( k = i + 1; k < m; k++ )
            {
              temp = temp - a[k+i*lda] * b[k+j*ldb];
            }
            if ( nounit )
            {
              temp = temp / a[i+i*lda];
            }
            b[i+j*ldb] = temp;
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
    if ( cblas_lsame ( transa, 'N' ) )
    {
      if ( upper )
      {
        for ( j = 0; j < n; j++ )
        {
          if ( alpha != 1.0 )
          {
            for ( i = 0; i < m; i++ )
            {
              b[i+j*ldb] = alpha * b[i+j*ldb];
            }
          }
          for ( k = 0; k < j; k++ )
          {
            if ( a[k+j*lda] != 0.0 )
            {
              for ( i = 0; i < m; i++ )
              {
                b[i+j*ldb] = b[i+j*ldb] - a[k+j*lda] * b[i+k*ldb];
              }
            }
          }
          if ( nounit )
          {
            temp = 1.0 / a[j+j*lda];
            for ( i = 0; i < m; i++ )
            {
              b[i+j*ldb] = temp * b[i+j*ldb];
            }
          }
        }
      }
      else
      {
        for ( j = n - 1; 0 <= j; j-- )
        {
          if ( alpha != 1.0 )
          {
            for ( i = 0; i < m; i++ )
            {
              b[i+j*ldb] = alpha * b[i+j*ldb];
            }
          }
          for ( k = j + 1; k < n; k++ )
          {
            if ( a[k+j*lda] != 0.0 )
            {
              for ( i = 0; i < m; i++ )
              {
                b[i+j*ldb] = b[i+j*ldb] - a[k+j*lda] * b[i+k*ldb];
              }
            }
          }
          if ( nounit )
          {
            temp = 1.0 / a[j+j*lda];
            for ( i = 0; i < m; i++ )
            {
              b[i+j*ldb] = temp * b[i+j*ldb];
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
      if ( upper )
      {
        for ( k = n - 1; 0 <= k; k-- )
        {
          if ( nounit )
          {
            temp = 1.0 / a[k+k*lda];
            for ( i = 0; i < m; i++ )
            {
              b[i+k*ldb] = temp * b[i+k*ldb];
            }
          }
          for ( j = 0; j < k; j++ )
          {
            if ( a[j+k*lda] != 0.0 )
            {
              temp = a[j+k*lda];
              for ( i = 0; i < m; i++ )
              {
                b[i+j*ldb] = b[i+j*ldb] - temp * b[i+k*ldb];
              }
            }
          }
          if ( alpha != 1.0 )
          {
            for ( i = 0; i < m; i++ )
            {
              b[i+k*ldb] = alpha * b[i+k*ldb];
            }
          }
        }
      }
      else
      {
        for ( k = 0; k < n; k++ )
        {
          if ( nounit )
          {
            temp = 1.0 / a[k+k*lda];
            for ( i = 0; i < m; i++ )
            {
              b[i+k*ldb] = temp * b[i+k*ldb];
            }
          }
          for ( j = k + 1; j < n; j++ )
          {
            if ( a[j+k*lda] != 0.0 )
            {
              temp = a[j+k*lda];
              for ( i = 0; i < m; i++ )
              {
                b[i+j*ldb] = b[i+j*ldb] - temp * b[i+k*ldb];
              }
            }
          }
          if ( alpha != 1.0 )
          {
            for ( i = 0; i < m; i++ )
            {
              b[i+k*ldb] = alpha * b[i+k*ldb];
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
void dsyrk(char uplo, char trans, int n, int k,
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
	else  if (lda < max(1, nrowa)) {
		info = 7;
	}
	else  if (ldc < max(1, n)) {
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
		/*
		*  Form  C := alpha*A*A**T + beta*C.
		*/
		if (upper) {
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
						temp = alpha*A[i * lda + l];
						for (i = 0; i < j; i++) {
							C[i * ldc + j] = C[i * ldc + j] + temp*A[i * lda + l];
						}
					}
				}
			}
		}
		else {

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
					if (A[i * lda + l] != ZERO) {
						temp = alpha*A[i * lda + l];
						for (i = j; i < n; i++) {
							C[i * ldc + j] = C[i * ldc + j] + temp*A[i * lda + l];
						}
					}
				}
			}
		}
	}
	else {
		/*
		*  Form  C := alpha*A**T*A + beta*C.
		*/
		if (upper) {
			for (j = 0; j < n; j++) {
				for (i = 0; i < j; i++) {
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
void cblas_dpotrf2(char uplo, int 	n, double *A, int lda, int info)
{
	double ONE = 1.0;
	double ZERO = 0.0;

	int i, j, nb;
	int jb;
	/*
	*     Test the input parameters
	*/
	info = 0;
	bool upper = cblas_lsame(uplo, 'U');
	if (!upper &  !cblas_lsame(uplo, 'L')) {
		info = -1;
	}
	else if (n < 0) {
		info = -2;
	}
	else if (lda < max(1, n)) {
		info = -4;
	}
	if (info != 0) {
		xerbla("DPOTRF2", -info);
		return;
	}

	/*
	*     Quick return if possible
	*/
	if (n != 0) {
		return;
	}

	/*
	*     Determine the block size for this environment.
	*/
	nb = ilaenv(1, "DPOTRF", uplo, n, -1, -1, -1);

	if ((nb <= 1) || (nb >= n)) {
		/*
		* Use unblocked code.
		*/
		cblas_dpotrf2(uplo, n, A, lda, info);
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
				jb = min(nb, n - j + 1);
				dsyrk('U', 'T', jb, j - 1, -ONE, &A[1 * lda + j], lda, ONE, &A[j * lda + j], lda);
				cblas_dpotrf2('U', jb, &A[j * lda + j], lda, info);
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
				jb = min(nb, n - j + 1);
				dsyrk('L', 'N', jb, j - 1, -ONE, &A[j * lda + 1], lda, ONE, &A[j * lda + j], lda);
				cblas_dpotrf2('L', jb, &A[j * lda + j], lda, info);
				if (info != 0) {
					goto MARK30;
				}
				if ((j + jb) <= n) {
					/*
					* Compute the current block column.
					*/
					cblas_dgemm('N', 'T', n - j - jb + 1, jb, j - 1, -ONE, &A[(j + jb)*lda + 1], lda, &A[j * lda + 1], lda, ONE, &A[(j + jb) * lda + j], lda);
					cblas_dtrsm('R', 'L', 'T', 'N', n - j - jb + 1, jb, ONE, &A[j * lda + j], lda, &A[(j + jb)*lda + j], lda);
				}
			}
		}
	}
	goto MARK40;
MARK30:
	info = info + j - 1;

MARK40:
	return;
}


