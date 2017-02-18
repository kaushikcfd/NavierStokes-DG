#include "DG_Element_2d.h"
#include "../Utilities/LobattoNodes.hpp"


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the constructor which initializes the member variables of the class. It also populates the arrays
 * X, Y as once the co-ordinates of the rectangle is known these arrays can be computed.
 *
 * @Param[in] _N The order of the polynomial which is to be used.
 * @Param[in] x1 The x-coord of the bottom left corner.
 * @Param[in] y1 The y-coord of the bottom left corner.
 * @Param[in] x2 The x-coord of the top right corner.
 * @Param[in] y2 The y-coord of the top right corner.
 */
/* ----------------------------------------------------------------------------*/
DG_Element_2d::DG_Element_2d(int _N, float x1, float y1, float x2, float y2) {
    N = _N;
    x_start = x1;
    y_start = y1;
    x_end   = x2;
    y_end   = y2;


    X = new float[(N+1)*(N+1)];
    Y = new float[(N+1)*(N+1)];

    float* nodes = new float[N+1]; // Alloting space for the lobatto nodes.
    
    lobattoNodes(nodes, N+1); // Found the lobatto nodes in the range -1 to +1.

    for(int k1=0;k1<=N;k1++)
    {
        for(int k2=0;k2<=N;k2++)
        {
            X[k1*(N+1)+k2] = 0.5*(x1 + x2) + 0.5*(x2 - x1)*nodes[k2];
            Y[k1*(N+1)+k2] = 0.5*(y1 + y2) + 0.5*(y2 - y1)*nodes[k1];
        }
    }

    delete[] nodes;

}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This functions creates space in order to take in one more variable on which operators are needed to be
 * applied.
 *
 * @Param v  This is the name of the variable which is to be added.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::addVariable_withoutBoundary(string v) {
    float * newVariable = new float[(N+1)*(N+1)]; /// Allocating the space for the new variable which is to be created.
    variable[v] = newVariable; /// Now assigning the same to the map.

    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This functions creates space in order to take in one more variable on which operators are needed to be
 * applied.
 *
 * @Param v  This is the name of the variable which is to be added.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::addVariable_withBoundary(string v) {
    float * newVariable = new float[(N+1)*(N+1)]; /// Allocating the space for the new variable which is to be created.
    variable[v] = newVariable; /// Now assigning the same to the map.
    
    // **b_top is used because it will store the address to the boundary element, So whenever the actual value in the float* of the variable is changed then this will also change automatically. The same holds for all the other following mentioned variables.
    
    float **b_top       = new float* [N+1]; 
    float **b_bottom    = new float* [N+1];
    float **b_left      = new float* [N+1];
    float **b_right     = new float* [N+1];

   
    float **n_top       = new float* [N+1];
    float **n_left      = new float* [N+1];
    float **n_right     = new float* [N+1];
    float **n_bottom    = new float* [N+1];

    for(int i=0; i<=N; i++){
        b_bottom[i] = &(variable[v][i]);    
        b_top[i]    = &(variable[v][N*(N+1)+i]);
        b_left[i]   = &(variable[v][i*(N+1)]);
        b_right[i]  = &(variable[v][i*(N+1)+N]);
    }


    boundaryTop[v]      = b_top;
    boundaryRight[v]    = b_right;
    boundaryBottom[v]   = b_bottom;
    boundaryLeft[v]     = b_left;


    neighboringTop[v]   = n_top;
    neighboringRight[v] = n_right;
    neighboringBottom[v]= n_bottom;
    neighboringLeft[v]  = n_left;

    boundaryVariables.push_back(v);

    return ;
}
