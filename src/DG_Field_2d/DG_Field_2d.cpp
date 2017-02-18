/**
 * @file DG_Field_2d.cpp
 * @Synopsis  This is the file for the class `Field`, which the data-structure for storing DG Fields
 * @author Kaushik Kulkarni
 * @version 1.0
 * @date 2017-02-18
 */

#include "DG_Field_2d.h"


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the constructor method which takes the following inputs. Once it gets the inputs, it starts
 * creating the Elements, which in turn starts allocating space to it, and also giving the inputs of domain information
 * to each and every element.
 *
 * @Param _nex The number of elements in the x-direction for the field
 * @Param _ney The number of elements in the y-direction for the field
 * @Param _N The order of interpolation of the polynomial which is to be used.
 * @Param _x1 The grid is a structured rectangular grid, and hence this corresponds to the x-coord of the bottom left
 * corner of the grid.
 * @Param _y1 This corresponds to the y-coordinate of the bottom left corner of the grid
 * @Param _x2 This corresponds to the x-coordinate of the top right corner of the grid
 * @Param _y2 This corresponds to the y- coordinate of the top right corner of the grid
 */
/* ----------------------------------------------------------------------------*/
DG_Field_2d::DG_Field_2d(int _nex, int _ney, int _N, float _x1, float _y1, float _x2, float _y2) {
    ne_x = _nex;
    ne_y = _ney;
    N = _N;
    x1 = _x1;
    x2 = _x2;
    y1 = _y1;
    y2 = _y2;

    elements.resize(ne_x);

    float x_curr,y_curr, dx = (x2-x1)/ne_x, dy = (x2-x1)/ne_y;
    x_curr = x1;
    y_curr = y1;

    for(int i=0; i<ne_x; i++){
        y_curr = y1;
        for(int j=0; j<ne_y; j++){
            elements[i].push_back( new DG_Element_2d(N, x_curr, y_curr, x_curr + dx, y_curr + dy) );
            y_curr += dy;
        }
        x_curr += dx;
    }

}

void DG_Field_2d::addVariable_withBounary(string v) {
    
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->addVariable_withBoundary(v); // Adding the variable for the (i, j) th element.
       }
   }
    return ;
}


void DG_Field_2d::addVariable_withoutBounary(string v) {
    
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->addVariable_withoutBoundary(v); // Adding the variable for the (i, j) th element.
       }
   }
    return ;
}
