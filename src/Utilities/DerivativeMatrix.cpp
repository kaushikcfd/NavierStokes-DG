#include <functional>
#include <cstring>
#include "../../include/Utilities/LagrangePolynomials.h"
#include "../../include/Utilities/PolyEval.h"
#include "../../include/Utilities/PolyDeriv.h"
#include "../../include/Utilities/LobattoIntegration.h"
#include "../../include/Utilities/MassMatrix.h"

using namespace std;

void derivativeMatrix(float *DerivativeMatrix,unsigned N)
{
    float Poly[N+1][N+1];
    float **poly, **deriv;
    poly   =   new float*[N+1];
    deriv  =   new float*[N+1];
    lagrangePolynomials(*Poly,N);
    unsigned i,j;

    for(i=0;i<=N;i++)
    {
        poly[i] =   new float[N+1];
        deriv[i]=   new float[N];
        memcpy(poly[i],Poly[i],(N+1)*sizeof(float));
        polyDeriv(poly[i],deriv[i],N+1);
    }

    function<float(float)> eval;
    for(i=0;i<=N;i++)
    {
        for(j=0;j<=N;j++)
        {
            eval = [&poly,&deriv,&i,&j,&N](float x){return (((polyEval(poly[i],N,x))*(polyEval(deriv[j],N-1,x))));};
            DerivativeMatrix[i*(N+1)+j] = lobattoIntegration(-1.0,1.0,N+1,eval);
        }
    }

    for(i=0;i<=N;i++)
    {
        delete[] poly[i];
        delete[] deriv[i];
    }

    delete[] poly;
    delete[] deriv;
    return ;
}

void twoDDerivativeMatrixX(float *DerivativeMatrix, unsigned N)
{
    float m[N+1][N+1];
    float d[N+1][N+1];
    derivativeMatrix(*d,N);
    massMatrix(*m,N);

    unsigned i1,i2,j1,j2;
    for(i1=0;i1<=N;i1++)
        for(j1=0;j1<=N;j1++)
            for(i2=0;i2<=N;i2++)
                for(j2=0;j2<=N;j2++)
                    DerivativeMatrix[(i1*(N+1)+j1)*(N+1)*(N+1)+i2*(N+1)+j2] = m[i1][i2]*d[j1][j2];

    return ;
}

void twoDDerivativeMatrixY(float *DerivativeMatrix, unsigned N)
{
    float m[N+1][N+1];
    float d[N+1][N+1];
    derivativeMatrix(*d,N);
    massMatrix(*m,N);

    unsigned i1,i2,j1,j2;
    for(i1=0;i1<=N;i1++)
        for(j1=0;j1<=N;j1++)
            for(i2=0;i2<=N;i2++)
                for(j2=0;j2<=N;j2++)
                    DerivativeMatrix[(i1*(N+1)+j1)*(N+1)*(N+1)+i2*(N+1)+j2] = d[i1][i2]*m[j1][j2];

    return ;
}