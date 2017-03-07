#include "../includes/DG_Field_2d/DG_Field_2d.h"
#include "../includes/Utilities/DerivativeMatrix.h"
#include "../includes/Utilities/MassMatrix.h"
#include "../includes/Utilities/FluxMatrix.h"
#include <iostream>
#include <cmath>

using namespace std;

double Q(double x, double y) {
    return (x+y*y);
}

double Q_yy_exact(double x, double y) {
    return (2.0);
}

int main() {
    DG_Field_2d* field;
    field = new DG_Field_2d(10, 10, 16, -1.0, -1.0, 1.0, 1.0);
    
    field->addVariable_withBounary("Q");
    field->addVariable_withBounary("Q_y");
    field->addVariable_withBounary("Q_yy");
    field->addVariable_withoutBounary("Q_yy_exact");
    
    field->initializeVariable("Q", Q);
    field->initializeVariable("Q_yy_exact", Q_yy_exact);
    
    field->delByDelY("Q"  , "Q_y" , "central");
    field->delByDelY("Q_y" , "Q_yy", "central");

    cout << field->l2Norm("Q_yy_exact", "Q_yy") << endl;

    //field->writeVTK("output.vtk");
}
