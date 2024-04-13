#include <iostream>
#include "cblas.h"
#include "mdata.h"


void test_dlaswp()
{
	std::cout << std::endl;
	std::cout << "    Matrix Class Test  " << std::endl;
	std::cout << " ------------------------------------------ " << std::endl;
	std::cout << "    dlaswp  " << std::endl;
	std::cout << std::endl;

	double* pData_1 = new double[4 * 4]
	{
		1, 12, 13, 14,
		21, 22, 23, 24,
		31, 32, 33, 34,
		41, 42, 43, 44,
	};

	char trans = 'N';
	int n = 4;
	int nrhs = 1;
	double* a = pData_1;
	int lda = 4;
	int* ipiv = new int[n] {0, 0, 0, 0};
	//double* b = pB;
	int ldb = 4;
	int info = 0;
	int k1 = 0;
	int k2 = 2;
	int incx = 1;
	ipiv[0] = 2;
	ipiv[1] = 2;
	dlaswp( n, pData_1, lda,  k1, k2, ipiv, incx);


	MData theM(4, 4);
	theM.setData(a);

	theM.show();

}


void dgetrf2_test()
{
	printf("      cblas Library Test Case    \n");
	printf(" ----------------------------------- \n");
	printf("   dgetrf2 computes an LU factorization of a general M-by-N matrix A  \n");
	printf("	using partial pivoting with row interchanges.   \n");
	printf(" \n");

	int m = 3;
	int n = 3;
	int nrhs = 1; 
	int lda = 3;
	int info = 0;

	double* pData_1 = new double[3 * 3]
	{ 
		2, 4, 2,
		4, 8, 3,
		3, 4, 2
	};

	//double* pData_1 = new double[4 * 4]
	//{
	//	4, 3, 2, 8,
	//	21, 22, 23, 24,
	//	31, 32, 33, 34,
	//	41, 42, 43, 44,
	//};

	int* ipiv = new int[n] {0, 0, 0, 0};

	MData theM(3, 3);
	theM.setData(pData_1);
	theM.show();

	m = 3;
	n = 3; 
	lda = 3;
	dgetrf2(m, n, pData_1, lda, ipiv, info);

	theM.show();

	//char TRANS = 'N';
	//int LDA;
	//int INCX = 1;
	//int INCY = 1;

	//double* A;
	//int m, n, k;
	//double alpha = 1.0;
	//double beta = 1.0;
	//int* IPIV;
	//int INFO = 0;

	//int i;

	//m = 3, k = 2, n = 3;
	//k = std::min(m, n);


	//printf(" Initializing data \n");

	//printf(" Allocating memory for matrices \n\n");
	//A = (double*)malloc(m * n * sizeof(double));
	//IPIV = (int*)malloc(k * sizeof(int));

	//if (A == NULL || IPIV == NULL) {
	//	printf("\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
	//	free(A);
	//	free(IPIV);
	//	return;
	//}

	//for (i = 0; i < (m * n); i++) {
	//	A[i] = i * 1.0;
	//}
	//A[0] = 10.0;
	//A[1] = -7.0;
	//A[2] = 0.0;
	//A[3] = -3.0;
	//A[4] = 2.0;
	//A[5] = 6.0;
	//A[6] = 5.0;
	//A[7] = -1.0;
	//A[8] = 5.0;

	//for (i = 0; i < k; i++) {
	//	IPIV[i] = i * 1.0;
	//}


	//printf(" A: \n");
	////showMatrix_d(A, m, n);

	//dgetrf2(m, n, A, m, IPIV, INFO);


	//printf(" After Decomposition, A: \n");
	////showMatrix_d(A, m, n);
}




void test_dtrsm()
{
	std::cout << std::endl;
	std::cout << "    Matrix Class Test  " << std::endl;
	std::cout << " ------------------------------------------ " << std::endl;
	std::cout << "    DTRSM solves A*X=alpha*B or X*A=alpha*B, for triangular A, rectangular B.  " << std::endl;
	std::cout << std::endl;

	double* pData_1 = new double[4 * 4]
	{
		1, 12, 13, 14,
		21, 1, 3, 4,
		31, 32, 1, 4,
		41, 42, 43, 1,
	};

	double* pData_2 = new double[4]{ 1, 10, 1, 10 };

	char side = 'L';
	char uplo = 'L';
	char transa = 'N';
	char diag = 'U';
	int m = 4;
	int n = 1;
	double alpha = 1.0;
	double* pA = pData_1;
	int lda = 4;
	double* pB = pData_2;
	int ldb = 1;

	dtrsm(side, uplo, transa, diag, m, n,
		alpha, pA, lda, pB, ldb);

	//MData theM(4, 4);
	//theM.setData(pA);

	//theM.show();

	MData theB(4, 1);
	theB.setData(pB);

	theB.show();
	/*****
	* Results:
	1.000
	- 11.000
	322.000
	- 13415.000
	-- PASS --
	*/

	pB[0] = 1;
	pB[1] = 10;
	pB[2] = 1;
	pB[3] = 10;

	uplo = 'U';
	dtrsm(side, uplo, transa, diag, m, n,
		alpha, pA, lda, pB, ldb);
	theB.show();
	/**
	Results:
	-676.000
	87.000
	-39.000
	10.000
	-- PASS --
	*/


}



void dtrmm_test()
{
	std::cout << std::endl;
	std::cout << "    Matrix Class Test  " << std::endl;
	std::cout << " ------------------------------------------ " << std::endl;
	std::cout << "    DTRMM performs B:=A*B or B:=B*A, A triangular, B rectangular.  " << std::endl;
	std::cout << std::endl;

	double* pData_1 = new double[4 * 4]
	{
		1, 12, 13, 14,
		21, 1, 3, 4,
		31, 32, 1, 4,
		41, 42, 43, 1,
	};


	double* pData_2 = new double[4 * 4]
	{
		1, 1, 1, 1,
		2, 1, 3, 4,
		3, 2, 1, 4,
		4, 4, 3, 1,
	};

	char side = 'L';
	char uplo = 'U';
	char transa = 'N';
	char diag = 'U';
	int m = 4;
	int n = 4;
	double alpha = 1.0;
	double* pA = pData_1;
	int lda = 4;
	double* pB = pData_2;
	int ldb = 4;

	MData theA(4, 4);
	theA.setData(pData_1); 
	std::cout << "A:  " << std::endl;
	theA.show();

	MData theB(4, 4);
	theB.setData(pB);
	std::cout << "B:  " << std::endl;
	theB.show();


	dtrmm(side, uplo, transa, diag, m, n,
		alpha, pA, lda, pB, ldb);


	std::cout << "After trmm:  " << std::endl;
	theB.show();


}


void test_dtrsv()
{
	std::cout << std::endl;
	std::cout << "    Matrix Class Test  " << std::endl;
	std::cout << " ------------------------------------------ " << std::endl;
	std::cout << "    DTRSV  solves one of the systems of equations A*x = b, or A**T*x = b,  " << std::endl;
	std::cout << std::endl;

	double* pData_1 = new double[4 * 4]
	{
		1,  2,  3, 1 ,
		 1, 1, 3, 4,
		 1,  2, 1, 4,
		 1,  2,  3, 1,
	};

	double* pData_2 = new double[4]{ 1, 10, 1, 10 };

	char side = 'L';
	char uplo = 'U';
	char trans = 'T';
	char diag = 'U';
	int m = 4;
	int n = 4;
	double alpha = 1.0;
	double* pA = pData_1;
	int lda = 4;
	double* x = pData_2;
	int incx = 1;

	MData theM(4, 4);
	theM.setData(pA);

	std::cout << "Matrix A:" << std::endl;
	theM.show();

	MData theB(4, 1);
	theB.setData(x);
	std::cout << "b:" << std::endl;
	theB.show();

	dtrsv(uplo, trans, diag, n, pA, lda, x, incx);
	//dtrsv(side, uplo, transa, diag, m, n,
	//	alpha, pA, lda, pB, ldb);

	std::cout << "Result x:" << std::endl;
	theB.show();

	/*****
	  Results:
 
	-- PASS --
	*/


}




void dgetrf_test()
{
	printf("      cblas Library Test Case    \n");
	printf(" ----------------------------------- \n");
	printf("   dgetrf computes an LU factorization of a general M-by-N matrix A  \n");
	printf("	using partial pivoting with row interchanges.   \n");
	printf(" \n");

	int m = 3;
	int n = 3;
	int nrhs = 1;
	int lda = 3;
	int info = 0;

	double* pData_1 = new double[3 * 3]
	{
		2, 4, 2,
		4, 8, 3,
		3, 4, 2
	};

	//double* pData_1 = new double[4 * 4]
	//{
	//	4, 3, 2, 8,
	//	21, 22, 23, 24,
	//	31, 32, 33, 34,
	//	41, 42, 43, 44,
	//};

	int* ipiv = new int[n] {0, 0, 0, 0};

	MData theM(3, 3);
	theM.setData(pData_1);
	theM.show();

	m = 3;
	n = 3;
	lda = 3;
	dgetrf(m, n, pData_1, lda, ipiv, info);

	theM.show();

	//char TRANS = 'N';
	//int LDA;
	//int INCX = 1;
	//int INCY = 1;

	//double* A;
	//int m, n, k;
	//double alpha = 1.0;
	//double beta = 1.0;
	//int* IPIV;
	//int INFO = 0;

	//int i;

	//m = 3, k = 2, n = 3;
	//k = std::min(m, n);


	//printf(" Initializing data \n");

	//printf(" Allocating memory for matrices \n\n");
	//A = (double*)malloc(m * n * sizeof(double));
	//IPIV = (int*)malloc(k * sizeof(int));

	//if (A == NULL || IPIV == NULL) {
	//	printf("\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
	//	free(A);
	//	free(IPIV);
	//	return;
	//}

	//for (i = 0; i < (m * n); i++) {
	//	A[i] = i * 1.0;
	//}
	//A[0] = 10.0;
	//A[1] = -7.0;
	//A[2] = 0.0;
	//A[3] = -3.0;
	//A[4] = 2.0;
	//A[5] = 6.0;
	//A[6] = 5.0;
	//A[7] = -1.0;
	//A[8] = 5.0;

	//for (i = 0; i < k; i++) {
	//	IPIV[i] = i * 1.0;
	//}


	//printf(" A: \n");
	////showMatrix_d(A, m, n);

	//dgetrf2(m, n, A, m, IPIV, INFO);


	//printf(" After Decomposition, A: \n");
	////showMatrix_d(A, m, n);
}


void dsyrk_test()
{
	printf("      cblas Library Test Case    \n");
	printf(" ----------------------------------- \n");
	printf("   DSYRK performs one of the symmetric rank k operations \n");
	printf("	  C : = alpha * A * A * *T + beta * C,   \n");
	printf("      C := alpha*A**T*A + beta*C,  \n");
	printf(" \n");
 
	char uplo = 'L';
	char trans = 'T';
	
	int n = 4;
	int k = 4;
	double alpha = 1.0;
	int lda = 4;

	double beta = 1.0; 
	int ldc = 4; 

	double* pData_1 = new double[4 * 4]
	{
		 4, 1, 3,  4,
		 1,  2,  3,  1,
		 3,  3,  3,  4,
		 4,  1 , 4, 4,
	};

	double* pC = new double[n*n] { 0, 1, 1, 1};


	dsyrk(uplo, trans, n, k, alpha, pData_1, lda, beta, pC, ldc);

	MData C(n, n);
	C.setData(pC);

	C.show();

}


void dlarfg_test()
{
	printf("      cblas Library Test Case    \n");
	printf(" ----------------------------------- \n");
	printf("   dlarfg  generates a real elementary reflector H of order n, such that \n"); 
	printf("	H* (alpha) = (beta), H** T* H = I.   \n");
	printf("		(x)(0)   \n");
	printf(" \n");

	int m = 4;
	int n = 4;
	int nrhs = 1;
	int lda = 3;
	int info = 0;
	double alpha = 1.0;
	double tau = 1.0;

	double* pData_1 = new double[3 * 3]
	{
		2, 4, 2,
		4, 8, 3,
		3, 4, 2
	};

	double* pX = new double[4 * 4]
	{
		4, 3, 2, 8,
		21, 22, 23, 24,
		31, 32, 33, 34,
		41, 42, 43, 44,
	};

	int* ipiv = new int[n] {0, 0, 0, 0};

	MData theM(n, n);
	theM.setData(pX);
	theM.show();
	 

	n = 4; 
	int incx = 4; 
	//alpha = pX[0];
	dlarfg(n, alpha, pX + incx, incx, tau);

	theM.show();
	std::cout << "alpha: " << alpha  << std::endl;
	std::cout << "tau: " << tau << std::endl;

}

