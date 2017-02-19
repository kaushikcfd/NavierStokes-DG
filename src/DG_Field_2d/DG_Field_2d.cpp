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


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function adds the variable to each and every element. In this the boundary points are not
 * specifically stored in each element.
 *
 * @Param v This is a string which defines the variable name.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::addVariable_withBounary(string v) {
    
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->addVariable_withBoundary(v); // Adding the variable for the (i, j) th element.
       }
   }
   variableNames.push_back(v);
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is similar to DG_Field_2d::addVariable_withBoundary. Just the boundary points are not specificaly
 * stored for this variable.
 *
 * @Param v This is the name of the variable which is to be added.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::addVariable_withoutBounary(string v) {
    
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->addVariable_withoutBoundary(v); // Adding the variable for the (i, j) th element.
       }
   }
   variableNames.push_back(v);
   return ;
}

void DG_Field_2d::initializeVariable(string v, function<float(float, float)> f) {
    

   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
            elements[i][j]->initializeVariable(v, f); // Initializing the corresponding element by passing the same parameters to it.
       }
   }

    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function in order to write the data in the form of VTK file.
 *
 * @Param fileName This is the string fileName with which the file is to be saved.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::writeVTK(string fileName){
    ofstream pFile;
    pFile.open(fileName);
    
    int i, j, k, k1, k2;

    // Printing the preamble for the .vtk file.
    pFile << "# vtk DataFile Version 3.0\nNavier Stokes DG\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    // The information of the number of points.
    pFile << "POINTS\t" << (N+1)*(N+1)*ne_x*ne_y << "\tfloat\n";

    // Writing the point co-ordinates.
    for ( j = 0; j < ne_y; j++ )
        for ( i = 0; i < ne_x; i++ )
            for( k = 0; k < (N+1)*(N+1); k++ )
                pFile << elements[i][j]->X[k] << "\t" << elements[i][j]->Y[k] <<"\t"<<0 <<endl;

    pFile << "\n\n";

    // Specifying the information about the CELLS.
    pFile << "CELLS\t" << (N*N*ne_x*ne_y) <<"\t" << 5*(N*N*ne_x*ne_y) << endl;

    // Writing information about the structure of the cells.
    for ( i = 0; i < ne_y; i++ ) {
        for ( j = 0; j < ne_x; j++ ) {
            for( k1 = 0; k1 < N; k1++ ) {
                for ( k2 = 0; k2 < N; k2++ ) {
                    k   =   (i*ne_x+j)*(N+1)*(N+1) +   k1*(N+1)    +   k2;
                    pFile << 4 << "\t" << k << "\t" << k+1 << "\t" << k+N+2 << "\t" << k+N+1 << endl;
                }
            }
        }
    }
    pFile << "\n\n";

    // Specifying the information about the CELL TYPES.
    pFile << "CELL_TYPES " << (N*N*ne_x*ne_y) << endl;

    // `9` is the CELL TYPE CODE for specifying that it is a quad.
    for ( i = 0; i < (N*N*ne_x*ne_y); i++)
        pFile << "9\n";
    pFile << "\n\n";

    // Specifying the information about the values of the scalars.
    
    pFile << "POINT_DATA\t"<< (N+1)*(N+1)*ne_x*ne_y <<"\n";
    
    int noOfVars = variableNames.size(); // Getting the number of variables
    float* currentVariable;

    for(k1=0; k1<noOfVars; k1++) {
        pFile << "SCALARS\t"<< variableNames[k1] <<"\tfloat\nLOOKUP_TABLE default\n";
        
        // Writing the value of the POINT_DATA, for the variable[variableNames[k1]] 
        for ( j = 0; j < ne_y; j++ ){
            for ( i = 0; i < ne_x; i++ ) {
                currentVariable = elements[i][j]->variable[variableNames[k1]];
                for( k = 0; k < (N+1)*(N+1); k++ ) {
                    pFile << currentVariable[k] << endl;
                }
            }
        }
    }
    pFile.close(); // Closing the file.
    return ;
}
