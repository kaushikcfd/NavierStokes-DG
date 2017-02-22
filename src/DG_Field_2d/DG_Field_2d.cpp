/**
 * @file DG_Field_2d.cpp
 * @Synopsis  This is the file for the class `Field`, which the data-structure for storing DG Fields
 * @author Kaushik Kulkarni
 * @version 1.0
 * @date 2017-02-18
 */

#include "../../include/DG_Field_2d/DG_Field_2d.h"
#include "../../include/DG_Element_2d/DG_Element_2d.h"

#include "../../include/Utilities/Inverse.h"
#include "../../include/Utilities/DerivativeMatrix.h"
#include "../../include/Utilities/MassMatrix.h"
#include "../../include/Utilities/FluxMatrix.h"



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
    
    /// Setting the grid, by setting the elements. The elements are set by providing their end points for the quads.
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
    } // All the elements have been initialized.

    /// Setting the interaction between the elements by passing their neighboring elements addresses to each of the
    //elements.
    
    // Setting the top elements of each of the elements.
    for(int i = 0; i < ne_x; i++)
        for(int j=0; j < (ne_y-1); j++)
            elements[i][j]->setNeighboringElement('T', elements[i][j+1]);
        

    // Setting the right elements of each of the elements.
    for(int i = 0; i < (ne_x-1); i++)
        for(int j=0; j < (ne_y); j++)
            elements[i][j]->setNeighboringElement('R', elements[i+1][j]);

    // Setting the bottom elements of each of the elements.
    for(int i = 0; i < ne_x; i++)
        for(int j=1; j < ne_y; j++)
            elements[i][j]->setNeighboringElement('B', elements[i][j-1]);
    
    // Setting the left elements of each of the elements.
    for(int i = 1; i < ne_x; i++)
        for(int j=0; j < ne_y; j++)
            elements[i][j]->setNeighboringElement('L', elements[i-1][j]);
    // All the neighboring elements have been set, except for the elements at the boundary.
    
    /// Computing and passing the mass matrices, derivative and flux matrices to each and every element.
    /// The main use of this step is to ensure that each and every elements is not computing the same matrices.
    
    /// Assigning spaces to those matrices.
    float *massMatrix = new float[(N+1)*(N+1)*(N+1)*(N+1)];
    float *derivativeMatrix_x = new float[(N+1)*(N+1)*(N+1)*(N+1)] ; /// This is the matrix which is laid out in 1-D, this would help us to find t    he $\frac{d}{dx}$ of any term. 
    float *derivativeMatrix_y = new float[(N+1)*(N+1)*(N+1)*(N+1)] ; /// This is the matrix which is laid out in 1-D, this would help us to find t    he $\frac{d}{dy}$ of any term.
    float* fluxMatrix_top = new float[(N+1)*(N+1)*(N+1)*(N+1)] ; /// This is the flux matrix for the top edge.
    float* fluxMatrix_right = new float[(N+1)*(N+1)*(N+1)*(N+1)] ; // The Flux Matrix for the right edge.
    float* fluxMatrix_bottom = new float[(N+1)*(N+1)*(N+1)*(N+1)] ; /// This would be the flux term for the the bottom edge.
    float* fluxMatrix_left = new float[(N+1)*(N+1)*(N+1)*(N+1)] ; /// The Flux matrix for the left edge.
    float* massInverse = new float[(N+1)*(N+1)*(N+1)*(N+1)] ;
    
    /// Calling functions to compute the entries of the matrix.
    
    twoDMassMatrix(massMatrix, N);
    inverse(massMatrix,massInverse,(N+1)*(N+1));
    twoDDerivativeMatrixX(derivativeMatrix_x, N);
    twoDDerivativeMatrixY(derivativeMatrix_y, N);
    
    twoDFluxMatrix3(fluxMatrix_top, N);
    twoDFluxMatrix2(fluxMatrix_right, N);
    twoDFluxMatrix4(fluxMatrix_left, N);
    twoDFluxMatrix1(fluxMatrix_bottom, N);

    /// Assigning the computed matrices to each and every element.
    for(int i = 0; i < ne_x; i++)
        for(int j=0; j < ne_y; j++){
            elements[i][j]->setMassMatrix(massMatrix);
            elements[i][j]->setInverseMassMatrix(massInverse);
            elements[i][j]->setderivateMatrix_x(derivativeMatrix_x);
            elements[i][j]->setderivateMatrix_y(derivativeMatrix_y);
            elements[i][j]->setTopFluxMatrix(fluxMatrix_top);
            elements[i][j]->setRightFluxMatrix(fluxMatrix_right);
            elements[i][j]->setLeftFluxMatrix(fluxMatrix_left);
            elements[i][j]->setBottomFluxMatrix(fluxMatrix_bottom);
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
   for (int i=0; i < ne_x; i++ ){
       for (int j=0; j<ne_y; j++) {
           elements[i][j]->setVariableNeighbors(v); // This is essential so that the addresses of the neighbors are stored in each and every element.
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

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to operate the partial derivative of the variable w.r.t. x.
 *
 * @Param v The variable which is to be differentiated
 * @Param vDash The variable in which the differentiated value is to be stored.
 * @Param fluxType The numerical flux type which is to be implemented while computing the derivative.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::delByDelX(string v, string vDash, string fluxType) {
    for(int i = 0; i < ne_x; i++ )
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->delByDelX(v, vDash, fluxType);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to operate the partial derivative of the variable w.r.t. y.
 *
 * @Param v The variable which is to be differentiated
 * @Param vDash The variable in which the differentiated value is to be stored.
 * @Param fluxType The numerical flux type which is to be implemented while computing the derivative.
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::delByDelY(string v, string vDash, string fluxType) {
    for(int i = 0; i < ne_x; i++ )
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->delByDelY(v, vDash, fluxType);

    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @brief      This is used to extern the function `saxpy` to the class DG_Field_2d. $\mathbf{y} \mapsto a\mathbf{x} + \mathbf{y}$
 *
 * @param[in]  a     The coefficient `a` of `x`
 * @param[in]  x     The column vector `x`
 * @param[in]  y     The column vector `y`
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::axpy(float a, string x, string y) {

    for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->axpy(a, x, y);
    
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @brief      This is used to extern the function `sscal` to the class DG_Field_2d. $\mathbf{x} \mapsto a\mathbf{x}$
 *
 * @param[in]  a     The coefficient `a` of `x`
 * @param[in]  x     The column vector `x`
 * @param[in]  y     The column vector `y`
 */
/* ----------------------------------------------------------------------------*/
void DG_Field_2d::scal(float a, string x) {
    for(int i = 0; i < ne_x; i++)
        for(int j = 0; j < ne_y; j++)
            elements[i][j]->scal(a, x);

    return ;
}