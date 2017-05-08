<pre>
void cblas_dcopy(const int N, const double *X, const int incX, double *Y,  
        const int incY)  
{  
  
    if (N <= 0)  
    {  
        return;  
    }  
    //  
    //        code for unequal increments or equal increments  
    //          not equal to 1  
    //  
    int i;  
    if ((incX != 1) || (incY != 1))  
    {  
        int ix = 1;  
        int iy = 1;  
        if (incX < 0)  
        {  
            ix = (1 - N) * incX + 1;  
        }  
        if (incY < 0)  
        {  
            iy = (1 - N) * incY + 1;  
        }  
        for ( i = 0; i < N; i++)  
        {  
            Y[iy] = X[ix];  
            ix = ix + incX;  
            iy = iy + incY;  
        }  
        return;  
    }  
  
    //  
    //        code for both increments equal to 1  
    //  
    //  
    //        clean-up loop  
    //  
    int m = (N % 7);  
    if (0 != m)  
    {  
        for ( i = 0; i < m; i++)  
        {  
            Y[i] = X[i];  
        }  
        if (N < 7)  
        {  
            return;  
        }  
    }  
    int mp1 = m ;  
    for ( i = mp1; i < N; i += 7)  
    {  
        Y[i] = X[i];  
        Y[i + 1] = X[i + 1];  
        Y[i + 2] = X[i + 2];  
        Y[i + 3] = X[i + 3];  
        Y[i + 4] = X[i + 4];  
        Y[i + 5] = X[i + 5];  
        Y[i + 6] = X[i + 6];  
    }  
    return;  
  
}  
</pre>
