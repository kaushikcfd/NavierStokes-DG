#include <lapacke.h>

#ifndef PolyDeriv_HPP
#define PolyDeriv_HPP

void polyDeriv(float *GivenPoly, float *DerivedPoly, unsigned deg)
{
	for(int i = 1;i<=deg;i++)
        DerivedPoly[i-1] = i*GivenPoly[i];
	return ;
}
#endif
