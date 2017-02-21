#ifndef PLOT_H
#define PLOT_H

#include <string>
using namespace std;

void plot(float *X, float *Y, float N, string GraphTitle = "", string LegendTitle = "", string fileName = "");

/**
  * A similar 1-D plot but with the y-axis ranges.
  */
void plot(float *X, float *Y, float N,float Ymin, float Ymax, string GraphTitle = "", string LegendTitle = "", string fileName = "");

/**
 * This is a plot function which takes three 2D arrays and makes the surface plots of them.
 * @param X    The X-Domain.
 * @param Y    The Y-Domain.
 * @param Z    The Z-Domain.
 */
void plot(float *X, float *Y, float *Z, unsigned Nx, unsigned Ny,float X1, float X2,float Y1, float Y2,float Z1, float Z2, string Name);

#endif
