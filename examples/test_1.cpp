#include <stdlib.h>
#include <stdio.h>
#include "cblas.h"

double cor( int n, double *pX, double *pY )
{
	double result = 0.0; 	
	if( n <= 1 )
	{
		printf("cor(): n = %d invalid", n);
		return result;
	}
	if( 0 == pX )
	{
		printf("cor(): pX invalid", n);
		return result;
	}
	if( 0 == pY )
	{
		printf("cor(): pY invalid", n);
		return result;
	}
	//
	double sum_X = cblas_dasum( n, pX, 1 );
	double sum_Y = cblas_dasum( n, pY, 1 );
	double mean_X = sum_X / n;
	double mean_Y = sum_Y / n;

	double sqrXSum = 0.0;
	double sqrYSum = 0.0;
	for(int i = 0; i < n; i++ )
	{
		sqrXSum += (pX[i] - mean_X)*(pX[i] - mean_X);
		sqrYSum += (pY[i] - mean_Y)*(pY[i] - mean_Y);
	}

	double stdDev_X = sqrt( sqrXSum / (n - 1) );
	double stdDev_Y = sqrt( sqrYSum / (n - 1) );

	double zSum = 0.0; 
	for(int i = 0; i < n; i++ )
	{
		zSum += ( (pX[i] - mean_X)/ stdDev_X )*( (pY[i] - mean_Y)/ stdDev_Y );
	}
	result = zSum / (n-1);

	return result;
}

void main()
{	  
  double X[] = {15, 18, 21, 24, 27}; 
  double Y[] = {25, 25, 27, 31, 32}; 

	double r = cor( 5, X, Y);
	printf("corefficient r = %f \n", r);
}
