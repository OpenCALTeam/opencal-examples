// Conway's game of Life transition function kernel

#include <OpenCAL-CL/calcl2D.h>


#define DEVICE_Q 0

__kernel void lifeTransitionFunction(__CALCL_MODEL_2D)
{
	calclThreadCheck2D();
	
	int i = calclGlobalRow();
	int j = calclGlobalColumns();
	
	if(i==0 && j==0) 	
		printf("sum = %f \n", calclGetSum2Di(DEVICE_Q));

 	CALint sizeOfX_ = calclGetNeighborhoodSize();

	int sum = 0, n;

/*	for (n=1; n<sizeOfX_; n++)
		sum += calclGetX2Di(MODEL_2D, DEVICE_Q, i, j, n);

	if ((sum == 3) || (sum == 2 && calclGet2Di(MODEL_2D, DEVICE_Q, i, j) == 1))
		calclSet2Di(MODEL_2D, DEVICE_Q, i, j, 1);
	else
		calclSet2Di(MODEL_2D, DEVICE_Q, i, j, 0);*/
}
