#include <lapacke.h>

#ifndef Zeros_HPP
#define Zeros_HPP

void zeros(float *x, unsigned n)
{
    unsigned i;
    for(i=0; i<n; i++)
        x[i]    =   0.0;
}

void zeros(float *A, unsigned m, unsigned n)
{
    unsigned i;
    for(i=0;i<m*n;i++)
            A[i] = 0.0;
}


#endif
