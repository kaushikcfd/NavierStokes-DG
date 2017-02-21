#include <lapacke.h>
#include <functional>
#include "../../include/Utilities/LobattoNodes.h"
#include "../../include/Utilities/LobattoWeights.h"
#include "../../include/Utilities/PolyEval.h"

using namespace std;


float lobattoIntegration(float start, float end, unsigned N, function<float(float)> f) {
    N++;
    float *Nodes,*Weights,*Values;
    Nodes   =   new float[N];
    Weights =   new float[N];
    Values  =   new float[N];
    unsigned i;
    float integral=0.0;

    if(start>=end)
	{
		printf("ERROR: Please look into it there is some error in the given inputs for lobattoIntegration()\n");
		return 0.0;
	}

	lobattoNodes(Nodes,N);
    lobattoWeights(Weights,N);

    for (i=0;i<N;i++)
	   Nodes[i] = 0.5*(start+end) +  (0.5*(start-end))*Nodes[i];///Made a shift from the computational space to the physical space.


	for(i=0;i<N;i++)
		Values[i] =   f(Nodes[i]);

	for(i=0;i<N;i++)
		integral += (Values[i]*Weights[i]);

	return (0.5*(end-start)*(integral));
}