#include <lapacke.h>
#include <cstring>
#include <cstdio>

#ifndef Inverse_HPP
#define Inverse_HPP

void inverse(float *A, float *Ainv,unsigned N)
{
    memcpy(Ainv,A,N*N*sizeof(float));
    int info;
    int ipiv[N];
    info =  LAPACKE_sgetrf(LAPACK_ROW_MAJOR,N,N,Ainv,N,ipiv);
    if(info==0)
        info =  LAPACKE_sgetri(LAPACK_ROW_MAJOR,N,Ainv,N,ipiv);
    if(info!=0)
        fprintf(stderr,"The inverse of the matrix was unsuccesful.\n");
}

#endif
