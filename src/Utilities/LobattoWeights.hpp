#include <lapacke.h>
#include <functional>

#include "LobattoNodes.hpp"
#include "PolyEval.hpp"
#include "PolyDeriv.hpp"
#include "LegendrePolynomial.hpp"

#ifndef LobattoWeights_HPP
#define LobattoWeights_HPP

void lobattoWeights(float *Weights, unsigned N)
{
    float *Poly, *Nodes;
    Poly    =   new float[N];
    Nodes   =   new float[N];
	legendrePolynomial(Poly,N-1);
    lobattoNodes(Nodes,N);

    function<float(float)> Eval;
    Weights[0]  = 2.0/((N)*(N-1));

    Eval = [&Poly,&N](float x){ return(polyEval(Poly,N-1,x));};

	for(int i=1;i<N-1;i++)
		Weights[i]    =   2/((N*(N-1))*(Eval(Nodes[i]))*(Eval(Nodes[i])));
	Weights[N-1] = (2.0/((N)*(N-1)));
}

#endif
