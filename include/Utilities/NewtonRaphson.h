#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include <functional>
using namespace std;

float newtonRaphson(function<float(float)> func ,float x);

#endif