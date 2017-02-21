#include "src/DG_Field_2d/DG_Field_2d.h"
#include <iostream>
#include <cmath>

using namespace std;

float initialPressure(float x, float y) {
    return (x + y*y*y);
}

int main() {
    DG_Field_2d* field;
    field = new DG_Field_2d(10, 10, 2, -1.0, -1.0, 1.0, 1.0);
    
    field->addVariable_withBounary("Q");
    field->addVariable_withBounary("Qdash");
    
    field->initializeVariable("Q", initialPressure);
    
    field->delBydelX("Q", "Qdash", "central");

    field->writeVTK("danda.vtk");

    
}
