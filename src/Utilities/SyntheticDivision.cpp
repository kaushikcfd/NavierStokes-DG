#include <cstring>

void syntheticDivision(float *GivenPoly, unsigned deg, float root)
{
    int i;
    float *Result;
    Result = new float [deg];


	Result[deg-1] = GivenPoly[deg];
	for(i = deg-2;i>=0;i--)
		Result[i] = Result[i+1]*root + GivenPoly[i+1];

    memcpy(GivenPoly,Result,deg*(sizeof(float)));
    GivenPoly[deg] = 0;

    delete [] Result;
    return ;
}