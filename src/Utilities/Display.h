#include <cstdio>


#ifndef DISPLAY_H
#define DISPLAY_H

void display(float *x, unsigned n)
{
    unsigned i;
    for(i=0; i<n; i++)
        printf("%6.2f\n",x[i]);
}

void display(float *A, unsigned m, unsigned n)
{
    unsigned i,j;
    for(i=0; i<m; i++)
    {
        for(j=0;j<n;j++)
            printf("%6.2f\t",A[i*n + j]);
        printf("\n");
    }
}

#endif