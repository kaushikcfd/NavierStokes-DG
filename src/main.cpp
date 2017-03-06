#include "../includes/DG_Field_2d/DG_Field_2d.h"
#include "../includes/Utilities/DerivativeMatrix.h"
#include "../includes/Utilities/MassMatrix.h"
#include "../includes/Utilities/FluxMatrix.h"
#include <iostream>
#include <cmath>

using namespace std;

float Q(float x, float y) {
    return (10*x + 30*y);
}

int main() {
    DG_Field_2d* field;
    field = new DG_Field_2d(10, 10, 16, -1.0, -1.0, 1.0, 1.0);
    
    field->addVariable_withBounary("Q");
    field->addVariable_withBounary("Q_x");
    field->addVariable_withBounary("Q_y");
    
    field->initializeVariable("Q", Q);
    
    field->delByDelX("Q"  , "Q_x" , "central");
    field->delByDelY("Q"  , "Q_y" , "central");

    field->writeVTK("output.vtk");
}
