#include "driver.h"
#include <iostream>
#include <cmath>

using namespace std;

float initialPressure(float x, float y) {
    return (exp(-8*(x*x+y*y)));
}

int main() {
    DG_Field_2d* P;
    P = new DG_Field_2d(10, 10, 4, -1.0, -1.0, 1.0, 1.0);
    P->addVariable_withoutBounary("Pressure");
    P->initializeVariable("Pressure", initialPressure);
    P->writeVTK("danda.vtk");

    
}