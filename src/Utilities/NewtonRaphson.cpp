#include <functional>
#include <cmath>
#define epsi 0.001
using namespace std;

float newtonRaphson(function<float(float)> func ,float x)
{
	float deriv;
	float x_curr= x;
	float x_pre = x;
	do
	{
		x_pre = x_curr;
		deriv = (func(x_pre+epsi)-func(x_pre-epsi))/(2*epsi);
		x_curr = x_pre - (func(x_pre))/(deriv);
	}while(fabs(x_curr - x_pre) >= 1e-6);
	return x_curr ;
}