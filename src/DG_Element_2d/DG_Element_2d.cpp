#include "../../include/DG_Element_2d/DG_Element_2d.h"

#include "../../include/Utilities/Zeros.h"
#include "../../include/Utilities/Inverse.h"

#include "../../include/Utilities/FluxMatrix.h"
#include "../../include/Utilities/MassMatrix.h"
#include "../../include/Utilities/DerivativeMatrix.h"
#include "../../include/Utilities/LobattoNodes.h"



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

    massMatrix = NULL;

    derivativeMatrix_x  =   NULL;
    derivativeMatrix_y  =   NULL; 
    fluxMatrix_top      =   NULL;
    fluxMatrix_right    =   NULL;
    fluxMatrix_bottom   =   NULL;
    fluxMatrix_left     =   NULL;
    inverseMassMatrix   =   NULL;
    topNeighbor         =   NULL;
    rightNeighbor       =   NULL;
    leftNeighbor        =   NULL;
    bottomNeighbor      =   NULL;

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
        n_bottom[i] =   b_bottom[i] = &(variable[v][i]);    
        n_top[i]    =   b_top[i]    = &(variable[v][N*(N+1)+i]);
        n_left[i]   =   b_left[i]   = &(variable[v][i*(N+1)]);
        n_right[i]  =   b_right[i]  = &(variable[v][i*(N+1)+N]);
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

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to set the value of a variable with the help of a function.
 *
 * @Param v This is the name of the variable whose value is to be  set.
 * @Param f This is the f(x, y) which is used to initialize the variable.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::initializeVariable(string v, function<float(float, float)> f) {
    for(int i=0; i<((N+1)*(N+1)); i++)
        variable[v][i] = f(X[i], Y[i]);
    
    return ;
}


void DG_Element_2d::setVariableNeighbors(string v) {
    int j;
    if(topNeighbor!=NULL) {
        for( j = 0 ; j <= N; j++) 
            neighboringTop[v][j] = topNeighbor->boundaryBottom[v][j];
    }
    if(rightNeighbor!=NULL) {
        for( j = 0 ; j <= N; j++) 
            neighboringRight[v][j] = rightNeighbor->boundaryLeft[v][j];
    }
    if(leftNeighbor!=NULL) {
        for( j = 0 ; j <= N; j++) 
            neighboringLeft[v][j] = leftNeighbor->boundaryRight[v][j];
    }
    if(bottomNeighbor!=NULL) {
        for( j = 0 ; j <= N; j++) 
            neighboringBottom[v][j] = bottomNeighbor->boundaryTop[v][j];
    }
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function sets the boundary variables for an element. This must be called whenever a new variable with
 * boundary is added. 
 *
 * @Param v The variable whose information about the boundaries is to be stored
 * @Param type The type of the neighbor whose information is to be looked upon.
 * @Param neighbor The pointer to the element 
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setNeighboringElement(char type, DG_Element_2d* neighbor) {
    string currentVariable;
    switch(type) {
        case 't' : // `t` or `T` for top 
        case 'T' :
            topNeighbor = neighbor;
            break;
        case 'r' : // `r` or `R` for right
        case 'R' :
            rightNeighbor = neighbor;
            break;
        case 'b' : // `b` or `B` for bottom
        case 'B' :
            bottomNeighbor = neighbor;
            break;
        case 'l' : // `l` or `L` for left
        case 'L' :
            leftNeighbor = neighbor;
            break;
        default:
            cout << "WARNING!. No such neighbor type " << type << endl;
    }
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to get the variable `v` differentiated partially w.r.t. `x` and then store it in the
 * variable `vDash`. The function also takes `fluxType` as an input which would describe the numerical scheme that
 * should be used in order to obtain the derivative.
 *
 * @Param v         Variable which is to be differentiated.
 * @Param vDash     Variable in which the derivative is to be stored.
 * @Param fluxType  The type of flux that is to be used. eg "central"
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::delByDelX(string v, string vDash, string fluxType) {
    float dy = (y_end - y_start);
    float dx = (x_end - x_start);

    if(fluxType == "central") {
        float* numericalFlux        =   new float[(N+1)*(N+1)]; /// Creating a temporary new variable.
        float* auxillaryVariable    =   new float[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        zeros(numericalFlux, (N+1)*(N+1));                                                       
        for(int i=0; i<=N; i++){
            numericalFlux[i*(N+1)+N]    = 0.5*( *boundaryRight[v][i]    + *neighboringRight[v][i] ) ;   
            numericalFlux[i*(N+1)]    = 0.5*( *boundaryLeft[v][i]     + *neighboringLeft[v][i] ) ;  
        }
        /// vDash = -0.5*dy*D*v
        cblas_sgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dy, derivativeMatrix_x, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

        /// Adding the numeical Flux terms as necessary.
        cblas_sgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dy, fluxMatrix_right,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
        cblas_sgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dy, fluxMatrix_left,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

        /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
        cblas_sgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);

        delete[] numericalFlux;
        delete[] auxillaryVariable;
    }
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This is the function to get the variable `v` differentiated partially w.r.t. `y` and then store it in the
 * variable `vDash`. The function also takes `fluxType` as an input which would describe the numerical scheme that
 * should be used in order to obtain the derivative.
 *
 * @Param v         Variable which is to be differentiated.
 * @Param vDash     Variable in which the derivative is to be stored.
 * @Param fluxType  The type of flux that is to be used. eg "central"
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::delByDelY(string v, string vDash, string fluxType) {
    float dy = (y_end - y_start);
    float dx = (x_end - x_start);

    if(fluxType == "central") {
        float* numericalFlux        =   new float[(N+1)*(N+1)]; /// Creating a temporary new variable.
        float* auxillaryVariable    =   new float[(N+1)*(N+1)]; /// Creating a temporary new variable, auxiallary variable
        zeros(numericalFlux, (N+1)*(N+1));                                                       
        for(int i=0; i<=N; i++){
            numericalFlux[N*(N+1)+i]    = 0.5*( *boundaryTop[v][i]    + *neighboringTop[v][i] ) ;   
            numericalFlux[i]            = 0.5*( *boundaryBottom[v][i]     + *neighboringBottom[v][i] ) ;  
        }
        /// vDash = -0.5*dy*D*v
        cblas_sgemv(CblasRowMajor, CblasTrans,   (N+1)*(N+1), (N+1)*(N+1), -0.5*dx, derivativeMatrix_y, (N+1)*(N+1), variable[v],   1, 0, auxillaryVariable, 1);

        /// Adding the numeical Flux terms as necessary.
        cblas_sgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  0.5*dx, fluxMatrix_top,   (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);
        cblas_sgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1),  -0.5*dx, fluxMatrix_bottom,    (N+1)*(N+1), numericalFlux, 1, 1, auxillaryVariable, 1);

        /// Multiplying my Mass Inverse, this is the final step in getting the derivative.
        cblas_sgemv(CblasRowMajor, CblasNoTrans, (N+1)*(N+1),(N+1)*(N+1), 4.0/(dx*dy), inverseMassMatrix,(N+1)*(N+1), auxillaryVariable,1,0,variable[vDash],1);

        delete[] numericalFlux;
        delete[] auxillaryVariable;
    }
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  This function takes Mass Matrix as an input
 *
 * @Param m This is the massMatix array(actually a matrix, but implemented as a 1-d array).
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setMassMatrix(float *m) {
    massMatrix = m;
    return ;
}

void DG_Element_2d::setInverseMassMatrix(float* im) {
    inverseMassMatrix = im;
    return ;
}

/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to set the x-derivative matrix.
 *
 * @Param d The array of the x-derivative matrix which is given as an input.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setderivateMatrix_x(float *d) {
    derivativeMatrix_x = d;
    return ;
}


/* ----------------------------------------------------------------------------*/
/**
 * @Synopsis  Function to set the y-derivative matrix.
 *
 * @Param d The array of the y-derivative matrix which is given as an input.
 */
/* ----------------------------------------------------------------------------*/
void DG_Element_2d::setderivateMatrix_y(float *d) {
    derivativeMatrix_y = d;
    return ;
}

/***THE FOLLOWING 4 FUNCTIONS HAVE THE SIMILAR JOB TO SET THE FLUX MATRICES.****/

void DG_Element_2d::setTopFluxMatrix(float* f) {
    fluxMatrix_top = f;
    return ;
}
void DG_Element_2d::setRightFluxMatrix(float* f){
    fluxMatrix_right = f;
    return ;
}

void DG_Element_2d::setLeftFluxMatrix(float* f){
    fluxMatrix_left = f;
    return ;
}

void DG_Element_2d::setBottomFluxMatrix(float* f){
    fluxMatrix_bottom = f;
    return ;
}