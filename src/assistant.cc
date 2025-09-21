/*
* assistant.cc
*
*  Created on: 2023-05-28
*      Author: Yang
*/
#include "cblas.h" 
#include "mdata.h"
#include <limits>
#include <string.h>

//#define min(x,y) (((x) < (y)) ? (x) : (y))
//#define max(x,y) (((x) > (y)) ? (x) : (y)) 

extern double cblas_dnrm2(const int N, const double *X, const int incX);
extern void cblas_dscal(const int N, const double DA, double *DX, const int INCX);
extern int cblas_dgemv(char TRANS, int M, int N, double ALPHA, double* pA, int LDA,
	double* X, int INCX, double BETA, double* Y, int INCY);
extern int cblas_dger(int M, int N, double ALPHA, double* X, int INCX, double* Y,
	int INCY, double* A, int LDA);


#if defined(DBL_MAX)
#define _MAX_DOUBLE_  DBL_MAX
#else
#define _MAX_DOUBLE_ 1.7976931348623158e+308
#endif

#if defined(DBL_MIN)
#define _MIN_DOUBLE_  DBL_MIN
#else
#define _MIN_DOUBLE_  2.22507385850720200e-308
#endif

 
double sign(double A, double B)
{
	if (B >= 0)
	{
		return fabs(A);
	}
	else {
		return -fabs(A);
	}
}

/*
CMACH is CHARACTER*1
Specifies the value to be returned by DLAMCH:
= 'E' or 'e',   DLAMCH := EPS
= 'S' or 's ,   DLAMCH := SFMIN
= 'B' or 'b',   DLAMCH := base
= 'P' or 'p',   DLAMCH := EPS*base
= 'N' or 'n',   DLAMCH := t
= 'R' or 'r',   DLAMCH := rnd
= 'M' or 'm',   DLAMCH := emin
= 'U' or 'u',   DLAMCH := rmin
= 'L' or 'l',   DLAMCH := emax
= 'O' or 'o',   DLAMCH := rmax
where
EPS   = relative machine precision
SFMIN = safe minimum, such that 1/SFMIN does not overflow
base  = base of the machine
prec  = EPS*base
t     = number of (base) digits in the mantissa
rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
emin  = minimum exponent before (gradual) underflow
rmin  = underflow threshold - base**(emin-1)
emax  = largest exponent before overflow
rmax  = overflow threshold  - (base**emax)*(1-EPS)
*/
double dlamch(const char* CMACH)
{
	double ONE = 1.0;
	double ZERO = 0.0;
	double RND, EPS, SFMIN, SMALL, RMACH;
	bool LSAME;
	// int digits, epsilon, huge, maxexponent, minexponent, radix, tiny;

	//Assume rounding, not chopping.Always.
	RND = ONE;
	if (ONE == RND) {
		EPS = std::numeric_limits<double>::epsilon() * 0.5; // EPSilon(ZERO) * 0.5;
	}
	else {
		EPS = std::numeric_limits<double>::epsilon(); // EPSilon(ZERO);
	}
	//
	if (cblas_lsame(*CMACH, 'E')) {
		RMACH = std::numeric_limits<double>::epsilon(); // EPS;
	}
	else  if (cblas_lsame(*CMACH, 'S')) {
		SFMIN = std::numeric_limits<double>::epsilon(); // tiny(ZERO);
		SMALL = ONE / _MAX_DOUBLE_; //ONE / huge(ZERO);
		if (SMALL >= SFMIN) {
			//
			//           Use SMALL plus a bit, to avoid the possibility of rounding
			//           causing overflow when computing  1/SFMIN.
			//
			SFMIN = SMALL*(ONE + EPS);
		}
		RMACH = SFMIN;
	}
	else if (cblas_lsame(*CMACH, 'B')) {
		RMACH = 10; // radix(ZERO);
	}
	else  if (cblas_lsame(*CMACH, 'P')) {
		RMACH = EPS * 10; // EPS * radix(ZERO);
	}
	else if (cblas_lsame(*CMACH, 'N')) {
		RMACH = std::numeric_limits<double>::digits; // digits(ZERO);
	}
	else if (cblas_lsame(*CMACH, 'R')) {
		RMACH = RND;
	}
	else  if (cblas_lsame(*CMACH, 'M')) {
		RMACH = std::numeric_limits<double>::min_exponent ; // minexpONEnt(ZERO);
	}
	else  if (cblas_lsame(*CMACH, 'U')) {
		RMACH = std::numeric_limits<double>::epsilon(); // tiny(ZERO);
	}
	else if (cblas_lsame(*CMACH, 'L')) {
		RMACH = std::numeric_limits<double>::max_exponent ; // maxexpONEnt(ZERO)
	}
	else  if (cblas_lsame(*CMACH, 'O')) {
		RMACH = _MAX_DOUBLE_; // huge(ZERO);
	}
	else {
		RMACH = ZERO;
	}
	//
	return RMACH;
}


/*
DLAMC3  is intended to force  A  and  B  to be stored prior to doing
the addition of  A  and  B ,  for use in situations where optimizers
might hold one of these in a register.
*/
double dlamc3(double A, double B)
{
	return A + B;
}

/***
@brief IDAMAX finds the index of the first element having maximum absolute value.

\param[in] N: int, number of elements in input vector(s)

\param[in] DX: double*, dimension ( 1 + ( N - 1 )*abs( INCX ) )

\param[in] INCX: int, storage spacing between elements of DX
*/
double cblas_idamax(int n, double* dx, int incx)
{
	double dmax;
	int i, ix;

	// Quick return
	int result = 0;
	if (n < 1 || incx <= 0)
	{
		return result;
	}
	result = 0;
	if (n == 1)
	{
		return result;
	}
	if (incx == 1) {
		//
		//        code for increment equal to 1
		//
		dmax = fabs(dx[0]);
		for (i = 1; i < n; i++) {
			if (fabs(dx[i]) > dmax) {
				result = i;
				dmax = fabs(dx[i]);
			}
		}
	}
	else {
		//
		//  code for increment not equal to 1
		//
		ix = 0;
		dmax = fabs(dx[0]);
		ix = ix + incx;
		for (i = 1; i < n; i++) {
			if (fabs(dx[ix]) > dmax) {
				result = i;
				dmax = fabs(dx[ix]);
			}
			ix = ix + incx;
		}
	}
	return result;
}



// DLARFG generates a real elementary reflector H of order n, such
// that
//
//       H * ( alpha ) = ( beta ),   H**T * H = I.
//           (   x   )   (   0  )
//
// where alpha and beta are scalars, and x is an (n-1)-element real
// vector. H is represented in the form
//
//       H = I - TAU * ( 1 ) * ( 1 v**T ) ,
//                     ( v )
//
// where TAU is a real scalar and v is a real (n-1)-element
// vector.
//
// If the elements of x are all zero, then TAU = 0 and H is taken to be
// the unit matrix.
//
// Otherwise  1 <= TAU <= 2.
// \param[in] N
// \verbatim
//          N is INTEGER
//          The order of the elementary reflector.
// \endverbatim
//
// \param[in,out] ALPHA
// \verbatim
//          ALPHA is DOUBLE PRECISION
//          On entry, the value alpha.
//          On exit, it is overwritten with the value beta.
// \endverbatim
//
// \param[in,out] X
// \verbatim
//          X is DOUBLE PRECISION array, dimension
//                         (1+(N-2)*abs(INCX))
//          On entry, the vector x.
//          On exit, it is overwritten with the vector v.
// \endverbatim
//
// \param[in] INCX
// \verbatim
//          INCX is INTEGER
//          The increment between elements of X. INCX > 0.
// \endverbatim
//
// \param[out] TAU
// \verbatim
//          TAU is DOUBLE PRECISION
//          The value TAU.
// \endverbatim
// 
//  Authors:
//  ========
// 
// \author Univ. of Tennessee
// \author Univ. of California Berkeley
// \author Univ. of Colorado Denver
// \author NAG Ltd.
void cblas_dlarfg(int N, double& ALPHA, double* X, int INCX, double& TAU)
{
	double beta, xnorm, rsafmn, safmin;
	int knt;
	double zero = 0.0; 
	double one = 1.0;
	if (N < 1) {
		TAU = zero;
		return;
	}
	//
	xnorm = cblas_dnrm2(N - 1, X, INCX);
	//
	if (xnorm == zero) {
		//
		//        H  =  I
		//
		TAU = zero;
	}
	else {
		//
		//        general case
		//
		beta = -sign(dlapy2(ALPHA, xnorm), ALPHA);
		safmin = dlamch("S") / dlamch("E");
		knt = 0;
		if (fabs(beta) < safmin) {
			//
			//           XNORM, BETA may be inaccurate; scale X and recompute them
			//
			rsafmn = one / safmin;
		ROUTE10:       //continue;
			knt = knt + 1;
			cblas_dscal(N - 1, rsafmn, X, INCX);
			beta = beta*rsafmn;
			ALPHA = ALPHA*rsafmn;
			if ((abs(beta) < safmin) & (knt < 20))
			{
				//         GO TO 10
				goto ROUTE10;
			}

			//
			//           New BETA is at most 1, at least SAFMIN
			//
			xnorm = cblas_dnrm2(N - 1, X, INCX);
			beta = -sign(dlapy2(ALPHA, xnorm), ALPHA);
		}
		TAU = (beta - ALPHA) / beta;
		cblas_dscal(N - 1, one / (ALPHA - beta), X, INCX);
		//
		//        If ALPHA is subnormal, it may lose relative accuracy
		//
		for (int j = 0; j < knt; j++)
		{
			beta = beta*safmin;
		}
		ALPHA = beta;
	}
	//
	return;
}

/***
DLARF applies a real elementary reflector H to a real m by n matrix
C, from either the left or the right. H is represented in the form

H = I - tau * v * v**T

where tau is a real scalar and v is a real vector.

If tau = 0, then H is taken to be the unit matrix.
*/
void cblas_dlarf(char side, int m, int n, int l, double* v, int incv, double tau, double* pC, int ldc, double* work)
{
	double one = 1.0;
	double zero = 0.0;
	int i;

	int lastv = 0;
	int lastc = 0;
	bool applyleft = false;
	applyleft = cblas_lsame(side, 'L');

	MData c(n, ldc);
	c.setData(pC);

	if (tau != zero)
	{
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
		do {
			lastv = lastv - 1;
			i = i - incv;
		} while ((lastv > 0) && (v[i] == zero));

		if (applyleft) {
			//Scan for the last non - zero column in C(1:lastv, : ).
			lastc = iladlc(lastv, n, pC, ldc);
		}
		else {
			//Scan for the last non - zero row in C(:, 1 : lastv).
			lastc = iladlr(m, lastv, pC, ldc);
		}
	}

	//Note that lastc.eq.0 renders the BLAS operations null; no special
	//case is needed at this level.

	if (applyleft) {
		/*
		* Form  H* C
		*/
		if (lastv > 0) {
			/*
			* w(1:lastc, 1) := C(1:lastv, 1 : lastc) * *T * v(1:lastv, 1)
			*/
			cblas_dgemv('T', lastv, lastc, one, pC, ldc, v, incv, zero, work, 1);
			/*
			* C(1:lastv, 1 : lastc) : = C(...) - v(1:lastv, 1) * w(1:lastc, 1) * *T
			*/
			cblas_dger(lastv, lastc, -tau, v, incv, work, 1, pC, ldc);
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
			cblas_dgemv('N', lastc, lastv, one, pC, ldc, v, incv, zero, work, 1);
			/*
			*C(1:lastc, 1 : lastv) : = C(...) - w(1:lastc, 1) * v(1:lastv, 1) * *T
			*/
			cblas_dger(lastc, lastv, -tau, work, 1, v, incv, pC, ldc);
		}
	}

}



/*****
DLASV2 computes the singular value decomposition of a 2-by-2
triangular matrix
[  F   G  ]
[  0   H  ].
On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
right singular vectors for abs(SSMAX), giving the decomposition

[ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
[-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
*****/
void dlasv2(double f, double g, double h, double& ssmin, double& ssmax,
	double& snl, double& csl, double& snr, double& csr)
{
	double zero = 0.0;
	double half = 0.5;
	double one = 1.0;
	double two = 2.0;
	double four = 4.0;
	bool gasmal, swap;
	int pmax;
	double a, clt, crt, d, fa, ft, ga, gt, ha, ht, l, m, mm, r, s, slt, srt, t, temp, tsign, tt;

	ft = f;
	fa = abs(ft);
	ht = h;
	ha = abs(h);
	//
	//       PMAX points to the maximum absolute element of matrix
	//       PMAX = 1 if F largest in absolute values
	//       PMAX = 2 if G largest in absolute values
	//       PMAX = 3 if H largest in absolute values
	//
	pmax = 1;
	swap = (ha > fa);
	if (swap) {
		pmax = 3;
		temp = ft;
		ft = ht;
		ht = temp;
		temp = fa;
		fa = ha;
		ha = temp;
		// Now FA >= HA
	}
	gt = g;
	ga = abs(gt);
	if (ga == zero) {
		//  Diagonal matrix
		ssmin = ha;
		ssmax = fa;
		clt = one;
		crt = one;
		slt = zero;
		srt = zero;
	}
	else 
	{
		gasmal = true;
		if (ga > fa) {
			pmax = 2;
			if ((fa / ga) < dlamch("EPS")) {
				// Case of very large GA
				gasmal = false;
				ssmax = ga;
				if (ha > one) {
					ssmin = fa / (ga / ha);
				}
				else {
					ssmin = (fa / ga)*ha;
				}
				clt = one;
				slt = ht / gt;
				srt = one;
				crt = ft / gt;
			}
		}
		if (gasmal) {
			// Normal case
			d = fa - ha;
			if (d == fa) {
				// Copes with infinite F or H
				l = one;
			}
			else {
				l = d / fa;
			}
			// Note that 0 < L < 1
			//
			m = gt / ft;
			// Note that abs(M) .le. 1/macheps
			//
			t = two - l;
			// Note that T .ge. 1
			//
			mm = m*m;
			tt = t*t;
			s = sqrt(tt + mm);
			// Note that 1 .le. S .le. 1 + 1/macheps
			//
			if (l == zero) {
				r = abs(m);
			}
			else {
				r = sqrt(l*l + mm);
			}
			//  Note that 0 .le. R .le. 1 + 1/macheps
			//
			a = half*(s + r);
			//  Note that 1 .le. A .le. 1 + abs(M)
			//
			ssmin = ha / a;
			ssmax = fa*a;
			if (mm == zero) {
				// Note that M is very tiny
				//
				if (l == zero) {
					t = sign(two, ft)*sign(one, gt);
				}
				else {
					t = gt / sign(d, ft) + m / t;
				}
			}
			else {
				t = (m / (s + t) + m / (r + l))*(one + a);
			}
			l = sqrt(t*t + four);
			crt = two / l;
			srt = t / l;
			clt = (crt + srt*m) / a;
			slt = (ht / ft)*srt / a;
		}
	}
	if (swap) {
		csl = srt;
		snl = crt;
		csr = slt;
		snr = clt;
	}
	else {
		csl = clt;
		snl = slt;
		csr = crt;
		snr = srt;
	}
	//
	//     Correct signs of SSMAX and SSMIN
	//
	if (pmax == 1)
		tsign = sign(one, csr)*sign(one, csl)*sign(one, f);
	if (pmax == 2)
		tsign = sign(one, snr)*sign(one, csl)*sign(one, g);
	if (pmax == 3)
		tsign = sign(one, snr)*sign(one, snl)*sign(one, h);
	ssmax = sign(ssmax, tsign);
	ssmin = sign(ssmin, tsign*sign(one, f)*sign(one, h));
	return;
}

//double MAX(const double first, const double second)
//{
//	return (first >= second) ? first : second;
//}

/********
DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
overflow and unnecessary underflow.
*********/
double dlapy2(double x, double y)
{
	double result;
	double w, xabs, yabs, z, hugeval;
	bool x_is_nan, y_is_nan;

	double zero = 0.0;
	double half = 0.5;
	double one = 1.0;

	x_is_nan = isnan(x);
	y_is_nan = isnan(y);
	if (x_is_nan) result = x;
	if (y_is_nan) result = y;
	hugeval = dlamch("Overflow");
	//
	if (!(x_is_nan || y_is_nan)) {
		xabs = fabs(x);
		yabs = fabs(y);
		w = std::max(xabs, yabs);
		z = std::min(xabs, yabs);
		if( (z == zero) || (w > hugeval) ){
			result = w;
		}
		else {
			result = w*sqrt(one + (z / w)*(z / w));
		}
	}

	return result;
}



//
// XERBLA  is an error handler for the LAPACK routines.
// It is called by an LAPACK routine if an input parameter has an
// invalid value.  A message is printed and execution stops.
//
// Installers may consider modifying the STOP statement in order to
// call system-specific exception-handling facilities.
//
//    SRNAME is CHARACTER*(*)
//          The name of the routine which called XERBLA.
//    INFO is INTEGER
//          The position of the invalid parameter in the parameter list
//          of the calling routine.
void xerbla(const char* SRNAME, int INFO)
{
	//	std::cout << "CBLAS: On entry to " << SRNAME << " parameter number " << INFO
	//	<< " had an illegal value" << std::endl;
	printf("CBLAS: On entry to %s parameter number %i had an illegal value. \n", SRNAME, INFO);
	return;
}

void cblas_xerbla(const char* SRNAME, int INFO)
{
	//	std::cout << "CBLAS: On entry to " << SRNAME << " parameter number " << INFO
	//	<< " had an illegal value" << std::endl;
	printf("CBLAS: On entry to %s parameter number %i had an illegal value.\n", SRNAME, INFO);
	return;
}



/****
@brief ILAENV returns problem-dependent parameters for the local
environment.  See ISPEC for a description of the parameters.

https://netlib.org/lapack/explore-html-3.6.1/db/deb/tstiee_8f_a453ec05c38be7fc786a3489c72cd4999.html
In this version, the problem-dependent parameters are contained in
the integer array IPARMS in the common block CLAENV and the value
with index ISPEC is copied to ILAENV.  This version of ILAENV is
to be used in conjunction with XLAENV in TESTING and TIMING.

@param[in] ISPEC: int,
Specifies the parameter to be returned as the value of
ILAENV.
= 1: the optimal blocksize; if this value is 1, an unblocked
algorithm will give the best performance.
= 2: the minimum block size for which the block routine
should be used; if the usable block size is less than
this value, an unblocked routine should be used.
= 3: the crossover point (in a block routine, for N less
than this value, an unblocked routine should be used)
= 4: the number of shifts, used in the nonsymmetric
eigenvalue routines
= 5: the minimum column dimension for blocking to be used;
rectangular blocks must have dimension at least k by m,
where k is given by ILAENV(2,...) and m by ILAENV(5,...)
= 6: the crossover point for the SVD (when reducing an m by n
matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
this value, a QR factorization is used first to reduce
the matrix to a triangular form.)
= 7: the number of processors
= 8: the crossover point for the multishift QR and QZ methods
for nonsymmetric eigenvalue problems.
= 9: maximum size of the subproblems at the bottom of the
computation tree in the divide-and-conquer algorithm
(used by xGELSD and xGESDD)
=10: ieee NaN arithmetic can be trusted not to trap
=11: infinity arithmetic can be trusted not to trap

@param[in] NAME    char*[]
           The name of the calling subroutine, in either upper case or
           lower case.

@param[in] OPTS    char*[]
          The character options to the subroutine NAME, concatenated
          into a single character string.  For example, UPLO = 'U',
          TRANS = 'T', and DIAG = 'N' for a triangular routine would
          be specified as OPTS = 'UTN'.

@param[in] N1      int
@param[in] N2      int
@param[in] N3      int
@param[in] N4      int

Problem dimensions for the subroutine NAME; these may not all
     be required.

@return 
      >= 0: the value of the parameter specified by ISPEC
      < 0:  if ILAENV = -k, the k-th argument had an illegal value. 
*/
int ilaenv(int ispec, const char* name, char OPTS, int N1, int N2, int N3, int N4)
{
	if ((ispec >= 1) & (ispec <= 5)) {
/*
*        Return a value from the common block.
*/ 
		char t[5];
		char geqr[] = "GEQR ";
		char gelq[] = "GELQ ";
		strncpy(t, name, 5);
		int ret = strcmp(t, geqr);
		if (0 == ret)
		{
			if (N3 == 2) {
				return 2;
			}
			else {
				return 1;
			} 
		}
		ret = strcmp(t, gelq);
		if (0 == ret)
		{
			if (N3 == 2) {
				return 2;
			}
			else {
				return 1;
			}
		}
		return ispec;
	}
	else if (ispec == 6) {


	}

	return 0;
}



//int MAX(int first, int second)
//{
//	return first ? second : first > second;
//}

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

	INTA = (unsigned char)CA;
	INTB = (unsigned char)CB;

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


/**
@brief
ILADLC scans a matrix for its last non-zero column.
*/
int iladlc(int m, int n, double* pA, int lda)
{
	double zero = 0.0;
	int i;
	int result;

	MData a(m, lda);
	a.setData(pA);

	// Quick test for the common case where one corner is non - zero.
	if (n == 0) {
		result = n;
	}
	else if ((a(0, n - 1) != zero) || (a(m - 1, n - 1) != zero)) {
		result = n;
	}
	else {
		// Now scan each column from the end, returning with the first non - zero.
		for (result = n - 1; result >= 0; result--) {
			for (i = 0; i < m; i++) {
				if (a(i, result) != zero) {
					return result;
				}
			}
		}
	}
	return result;
}

/***
@brief
ILADLR scans a matrix for its last non-zero row.
@result:  -1 -- all the entries of the matrix are ZERO.
  0 -- only the first row is non-zero row. 
*/
int iladlr(int m, int n, double* pA, int lda)
{
	double zero = 0.0;
	int i, j;
	int result;

	MData a(m, n, pA, lda);

	//Quick test for the common case where one corner is non - zero.
	if (m == 0) {
		result = m;
	}
	else if ((a(m, 0) != zero) || (a(m, n) != zero)) {
		result = m;
	}
	else {
		// Scan up each column tracking the last zero row seen.
		result = -1;
		for (j = 0; j < n; j++) {
			i = m - 1;
			//do {
			//	i = i - 1;
			//} while ((a(std::max(i, 0), j) == zero) && (i >= 0));
			while ((a(std::max(i, 0), j) == zero) && (i >= 0))
			{
				i = i - 1;
			} 
			result = std::max(result, i);
		}
	}
	return result;
}


// Assistant Methods
/***
  @brief Show vector 
  @param[in] data: double*, data pointer 
  @param[in] m: int, vector dimension 
*/
void showVector_d(double* data, int m)
{
	if (NULL == data)
	{
		return;
	}
	for (int i = 0; i < m; i++) {
		printf("%12.5G", data[i]);
		if (0 != i)
		{
			if (0 == (i % 6)) {
				printf("\n");
			}
		}
	}
	printf("\n");
}


/***
  @brief Show matrix
  @param[in] data: double*, data pointer
  @param[in] m: int, matrix column number 
  @param[in] k: int, matrix row number
*/
void showMatrix_d(double* data, int m, int k)
{
	if (NULL == data)
	{
		return;
	}
	for (int i = 0; i<std::min(m, 6); i++) {
		for (int j = 0; j< std::min(k, 6); j++) {
			printf("%12.5G", data[j + i*k]);
		}
		printf("\n");
	}
}


/***
@brief
ILADLR scans a matrix for its last non-zero row.
@result:  0 -- all the entries of the matrix are ZERO.
  1 -- only the first row is non-zero row.
*/
int cblas_iladlr(int m, int n, double* pA, int lda)
{
	double zero = 0.0;
	int i, j;
	int result = -1;

	MData a(m, n, pA, lda);

	//Quick test for the common case where one corner is non - zero.
	if (m < 0) {
		result = -1;
	}
	else if ((a(m - 1, 0) != zero) || (a(m - 1, n - 1) != zero)) {
		result = m - 1;
	}
	else {
		// Scan up each column tracking the last zero row seen.
		result = -1;
		for (j = 0; j < n; j++) {
			i = m - 1;
			//do {
			//	i = i - 1;
			//} while ((a(std::max(i, 0), j) == zero) && (i >= 0));
			while ((a(std::max(i, 0), j) == zero) && (i >= 0))
			{
				i = i - 1;
			}
			result = std::max(result, i);
		}
	}
	return result + 1;
}


