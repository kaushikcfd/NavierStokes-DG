#include <lapacke.h>
#include <algorithm>
#include <functional>

#include "LegendrePolynomial.h"
#include "PolyEval.h"
#include "PolyDeriv.h"
#include "SyntheticDivision.h"
#include "NewtonRaphson.h"
using namespace std;

#ifndef LOBATTONODES_H
#define LOBATTONODES_H

void lobattoNodes(float *Nodes, unsigned N)
{
    float *Poly,*DerivedPoly;
    Poly        =   new float[N];
    DerivedPoly =   new float[N-1];
    float root;
	float InitialGuess=-1.0;
	unsigned deg=N-2;///To store the degree at each state of the polynomial.
    legendrePolynomial(Poly,N-1);
    polyDeriv(Poly,DerivedPoly,N-1);

    function<float(float)> eval;
	Nodes[0]   =   -1;
	for(int i=1;i<=N-2;i++)
	{
		eval = [&DerivedPoly,&deg](float x){ return(polyEval(DerivedPoly,deg,x));};
		root = newtonRaphson(eval,InitialGuess);
		syntheticDivision(DerivedPoly,deg,root);
        --deg;
		Nodes[i]  =   root;
	}
	Nodes[N-1] =   1;
	sort(Nodes,Nodes+N);

    delete[] Poly;
    delete[] DerivedPoly;
}

#endif