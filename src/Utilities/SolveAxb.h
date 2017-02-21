#include <lapacke.h>
#include <cstring>
#include <iostream>
using namespace std;

#ifndef SOLVEAXB_H
#define SOLVEAXB_H

void solveAxb(float *A, float *x, float *b, unsigned N)
{
    float *B    = new float[N*N];
    float *c    = new float[N];
    memcpy(B,A,N*N*sizeof(float));
    memcpy(c,b,N*sizeof(float));
    int ipiv[N];
    int info;

    info    =   LAPACKE_sgesv(LAPACK_ROW_MAJOR,N,1,B,N,ipiv,c,1);
    if(info!=0)
        cerr << "The Linear solve `Ax=b` was not succesful.\n";
   
    memcpy(x,c,N*sizeof(float));

    delete[] B;
    delete[] c;
    return ;
}

#endif