#ifndef LOBATTOINTEGRATION_H
#define LOBATTOINTEGRATION_H

#include <functional>

using namespace std;

float lobattoIntegration(float start, float end, unsigned N, function<float(float)> f);

#endif