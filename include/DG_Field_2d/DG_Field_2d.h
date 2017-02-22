/**
  * @file DG_Field_2d.h
  * @Synopsis  This is the file for the class `Field`, which the data-structure for storing DG Fields, This is the
  * header file and it only stores the declarations to all the functions and the member variables of the class.
  * @author Kaushik Kulkarni
  * @version 1.0
  * @date 2017-02-18
  */

#ifndef DG_FIELD_2D_H_
#define DG_FIELD_2D_H_

#include <vector>
#include <functional>
#include <fstream>

using namespace std; 

class DG_Element_2d; 

class DG_Field_2d {
private:
    int N; /// The order of the polynomial of interpolation.
    int ne_x, ne_y; /// Since this is a structured grid we can define the grid with the help of number of elements in the x-direction and in the y-direction.
    float x1, y1, x2, y2;


public:
   vector< vector<DG_Element_2d*> > elements;
   vector<string> variableNames; // This is stores all the variables which have been added to the field.
   vector<string> variablesWithBoundaryInfo; // This stores all the variables whose boundary info. is also known.

    DG_Field_2d(int _nex, int _ney, int _N, float _x1, float _y1, float _x2, float _y2);
    void addVariable_withBounary(string v);
    void addVariable_withoutBounary(string v);
    void initializeVariable(string v, function<float(float, float)>);
    void writeVTK(string fileName);
    
    // Operators on the field.
    void delByDelX(string v, string vDash, string fluxType);
    void delByDelY(string v, string vDash, string fluxType);

    // Functions to apply linear operations on the variables.
    void axpy(float a, string x, string y);
    void scal(float a, string x);
};

#endif
