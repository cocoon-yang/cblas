#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
#include "lapacke.h"
#include "cblas.h"
#include "foo.h"

double precision = 0.0001;

//
double func(const double* x, int n) //函数
{
    if(0 == x )
    {
        printf("func(): Invalid parameters\n");
        return 0.0;
    }
    return  x[0]*x[0] + x[1]*x[1] + 1.0 * x[1];
}

// Gradient of the objective function
void dfunc(const double* x, double *df, int n) 
{
    if(0 == x )
    {
        printf("dfunc(): Invalid parameters\n");
        return;
    }  
    if(0 == df )
    {
        printf("dfunc(): Invalid parameters\n");
        return;
    }  

    df[0] = 2.0 * x[0];    
    df[1] = 2.0 * x[1] + 1.0 ;  

    double normal = cblas_ddot( n, df, 1, df, 1); 

    if( fabs(normal) > precision )
    {
        df[0] = df[0] / normal; 
        df[1] = df[1] / normal;    
    } 
     return;
}

void updateX( double* x, const double *dir, double step_len, int n )
{
    int i = 0;
    if(0 == x )
    {
        printf("updateX(): Invalid parameters -- x\n");
        return;
    }  
    if(0 == dir )
    {
        printf("updateX(): Invalid parameters -- dir\n");
        return;
    }  

    for( i = 0; i < n; i++ )
    {
        x[i] += step_len * dir[i];
    }
    return;
}

//
//  paremeters:
//    itera_state -- int, 
//                    0:  iteration complete, 
//                 true:  iteration continue.
//          
int iteration( double* x, const double *df, double* direction, double step_len, int var_num, int *itera_state )
{
    if(0 == x )
    {
        printf("iteration(): Invalid parameters -- x\n");
        return 0;
    }  
    if(0 == df )
    {
        printf("iteration(): Invalid parameters -- dir\n");
        return 0;
    } 
    if(0 == direction )
    {
        printf("iteration(): Invalid parameters -- dir\n");
        return 0;
    } 
    if(0 == itera_state )
    {
        printf("iteration(): Invalid parameters -- itera_state\n");
        return 0;
    }  

    *itera_state = 1;
    //
    // direction[i] = df[i]
    // cblas_dcopy( var_num, df, 1, direction, 1);
    for (int i = 0; i < var_num; i++) {
        direction[i] = -df[i] ;
    }    

    // Logging  -- BEGIN -- 
    printf("iteration(): direction = \n");   
    for( int i = 0; i < var_num; i++ )
    {
        printf("%f \t", direction[i] );
    }
    printf("\n"); 
    // Logging  -- END -- 

    updateX( x, direction, step_len, var_num );

    // Logging  -- BEGIN -- 
    printf("iteration(): After updating, x = \n");   
    for( int i = 0; i < var_num; i++ )
    {
        printf("%f \t", x[i] );
    }
    printf("\n"); 
    // Logging  -- END -- 

    return 0;
}


int main()
{
    double ERROR = 0.0001;
    int itera_num = 20;
    int itera_state = 0;
    int var_num = 2;
    double step_len = 0.7;
    double obj = 0.0;

    double *x = (double*) malloc( var_num * sizeof(double));
    if(x == NULL)
    {
        printf("Allocing memory failed.\n");
        goto END;
    }

    double *df = (double*) malloc( var_num * sizeof(double));
    if(df == NULL)
    {
        printf("Allocing memory failed.\n");
        goto END;
    }

    double *direct = (double*) malloc( var_num * sizeof(double));
    if(direct == NULL)
    {
        printf("Allocing memory failed.\n");
        goto END;
    }

    itera_state = 1;

    double x0 = 1.1;
    x[0] = x0;    
    x[1] = x0 + 3.0;

    double obj_old = obj = func( x, var_num );

    for( int i = 0; i < itera_num; i++ )
    {
        // Logging 
        printf("\n" ); 
        printf("Iteration %d \n", i);

        printf("Old objective = %f \n", obj_old ); 
        printf("Before updating, x = \n");   
        for( int i = 0; i < var_num; i++ )
        {
            printf("%f \t", x[i] );
        }
        printf("\n"); 

        // Calculating gradient of the objective funxtion
        dfunc( x, df, var_num) ;
    
        printf("df = \n");   
        for( int i = 0; i < var_num; i++ )
        {
            printf("%f \t", df[i] );
        }
        printf("\n"); 

        iteration( x, df, direct,  step_len, var_num, &itera_state );

        if( 0 == itera_state )
        {
            printf("Optimization Iteration Teiminated: 0 == itera_state \n");
            break;
        }

        obj = func( x, var_num );
        printf("Objective = %f \n", obj ); 

        if( fabs(obj - obj_old) < ERROR )
        {
            printf("Optimization Iteration Teiminated: dealt_Objective < ERROR \n");
            break;
        }

        obj_old = obj;
    }


END:
    free(direct);
    direct = 0;
    free(x);
    x = 0;
    free(df);
    df = 0;
    return 0;
}
