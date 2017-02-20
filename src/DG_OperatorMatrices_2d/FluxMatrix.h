#include <lapacke.h>
#include "../Utilities/LagrangePolynomials.hpp"
#include "../Utilities/Zeros.hpp"
#include "MassMatrix.h"

#ifndef FluxMatrix_HPP
#define FluxMatrix_HPP

void fluxMatrix(float *FluxMatrix, unsigned N)
{
    float Poly[N+1][N+1];
    float **poly;
    poly   =   new float*[N+1];
    lagrangePolynomials(*Poly,N);
    unsigned i,j;

    for(i=0;i<=N;i++)
    {
        poly[i] =   new float[N+1];
        memcpy(poly[i],Poly[i],(N+1)*sizeof(float));
    }

    function<float(float)> eval;
    for(i=0;i<=N;i++)
    {
        for(j=0;j<=N;j++)
        {
            eval = [&poly,&i,&j,&N](float x){return (((polyEval(poly[i],N,x))*(polyEval(poly[j],N,x))));};
            FluxMatrix[i*(N+1)+j] = eval(1)-eval(-1);
        }
    }

    for(i=0;i<=N;i++)
        delete[] poly[i];

    delete[] poly;
    return ;
}

void twoDFluxMatrix1(float *Flux1, unsigned N)
{
    unsigned i1,i2,j1,j2;
    float M[N+1][N+1];
    massMatrix(*M,N);
    zeros(Flux1,(N+1)*(N+1),(N+1)*(N+1));

    i1=i2=0;

    for(j1=0;j1<=N;j1++)
            for(j2=0;j2<=N;j2++)
                Flux1[(i1*(N+1)+j1)*(N+1)*(N+1)+i2*(N+1)+j2] = M[j1][j2];
}

void twoDFluxMatrix3(float *Flux3, unsigned N)
{
    unsigned i1,i2,j1,j2;
    float M[N+1][N+1];
    massMatrix(*M,N);
    zeros(Flux3,(N+1)*(N+1),(N+1)*(N+1));

    i1=i2=N;

    for(j1=0;j1<=N;j1++)
            for(j2=0;j2<=N;j2++)
                Flux3[(i1*(N+1)+j1)*(N+1)*(N+1)+i2*(N+1)+j2] = M[j1][j2];
}

void twoDFluxMatrix2(float *Flux2, unsigned N)
{
    unsigned i1,i2,j1,j2;
    float M[N+1][N+1];
    massMatrix(*M,N);
    zeros(Flux2,(N+1)*(N+1),(N+1)*(N+1));

    j1=j2=N;
    for(i1=0;i1<=N;i1++)
        for(i2=0;i2<=N;i2++)
                Flux2[(i1*(N+1)+j1)*(N+1)*(N+1)+i2*(N+1)+j2] = M[i1][i2];
}

void twoDFluxMatrix4(float *Flux4, unsigned N)
{
    unsigned i1,i2,j1,j2;
    float M[N+1][N+1];
    massMatrix(*M,N);
    zeros(Flux4,(N+1)*(N+1),(N+1)*(N+1));

    j1=j2=0;
    for(i1=0;i1<=N;i1++)
        for(i2=0;i2<=N;i2++)
                Flux4[(i1*(N+1)+j1)*(N+1)*(N+1)+i2*(N+1)+j2] = M[i1][i2];
}

#endif
