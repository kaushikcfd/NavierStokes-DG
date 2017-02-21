#ifndef POLYDERIV_H
#define POLYDERIV_H

void polyDeriv(float *GivenPoly, float *DerivedPoly, unsigned deg)
{
	for(int i = 1;i<=deg;i++)
        DerivedPoly[i-1] = i*GivenPoly[i];
	return ;
}
#endif
