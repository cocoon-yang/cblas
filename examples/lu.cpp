#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
#include "lapacke.h"
#include "cblas.h"
#include "foo.h"


void largeTest( int num, double val )
{
    int m = num;
    int n = num;
    int LDA = m; 
    double *pData = 0, *b = 0, *x =0;
    lapack_int* ipiv = 0;

    pData = (double*)malloc( m * n * sizeof(double));
    if(0 == pData)
    {
        printf("Allocing memory failed. \n");
        goto FREE_MEMORY;
    }
    b = (double*)malloc( m * sizeof(double));
    if(0 == b)
    {
        printf("Allocing memory failed. \n");
        goto FREE_MEMORY;
    }
    x = (double*)malloc( m * sizeof(double));  
    if(0 == x)
    {
        printf("Allocing memory failed. \n");
        goto FREE_MEMORY;
    }

    for( int i = 0; i < m*n; i++ )
    {
        pData[i] = 0.0;
    }

    for( int i = 0; i < m; i++ )
    {
        pData[i*LDA + i] = i + 1;
    }

    for( int i = 0; i < m - 1; i++ )
    {
        pData[i*LDA + i + 1] = i + 1 + 0.5;
    }

    for( int i = 1; i < m; i++ )
    {
        pData[i*LDA + i - 1] = i + 1 - 0.5;
    }   

    for( int i = 0; i < m; i++ )
    {
        b[i] = i;
        x[i] = 0.0 ;
    } 

    pData[(m-2)*LDA + (m-2) + 1] = val ;
    pData[(m-1)*LDA + (m-1) - 1] = val;

    ipiv = (lapack_int*) malloc( sizeof(lapack_int) * m );
    if(0 == ipiv)
    {
        printf("Allocing memory failed. \n");
        goto FREE_MEMORY;
    }

    // Calling dgetrf of the LAPACK library
    lapack_int info = LAPACKE_dgetrf( LAPACK_ROW_MAJOR, m, n, pData, LDA, ipiv );
    
    // x[] = b[];
    cblas_dcopy( n, b, 1, x, 1);
    
    for (int i = 0; i < n; i++) {
        int piv = ipiv[i] - 1;
        double tmp = x[i];
        x[i] = x[piv]  ;// b[ipiv[i]]; 
        x[piv] = tmp;

        for (int k = 0; k < i; k++){
            x[i] -= pData[i*LDA + k] * x[k]; 
        }
    }
    
    for (int i = n - 1; i >= 0; i--) {
        for (int k = i + 1; k < n; k++)
            x[i] -= pData[i*LDA + k] * x[k];
        x[i] = x[i] / pData[i*LDA + i];
    }
    
    // Showing the results
    printf("m = %d,  value = %f \n", m, val );
    printf("x = \t");   
    for( int i = 0; i < m; i++ )
    {
        printf("%f \t", x[i] );
    }
    printf("\n"); 

FREE_MEMORY:
    free(ipiv);
    ipiv = 0;
    free(x);
    x = 0;
    free(b);
    b = 0;    
    free(pData);
    pData = 0;
    
    return;
}



/* Main program */
int main( int argc, char* argv[] )
{
 
   int num = 10;
   double val = 9.0;
   double diff = 0.2;

   if( argc > 1 )
   {
       num  = atoi(argv[1]);
   }
   if( argc > 2 )
   {
       val  = atof(argv[2]);
   }
   if( argc > 3 )
   {
       diff  = atof(argv[3]);
   }
   num = num > 0 ? num : 10;
   
   largeTest( num, val - diff ); 

   return 1;
}  
