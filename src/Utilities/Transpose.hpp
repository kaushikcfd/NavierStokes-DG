#include <lapacke.h>
#include <cstring>

#ifndef Transpose_HPP
#define Transpose_HPP

void transpose(float *A, unsigned N)
{
    unsigned i,j;
    float B[N][N];
    for(i=0;i<N;i++)
        for(j=0;j<N;j++)
            B[j][i] =   A[i*N+j];

    memcpy(A,&B[0][0],N*N*sizeof(float));
    return ;
}

#endif
