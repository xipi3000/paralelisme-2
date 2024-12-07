#include <math.h>
#include <float.h>
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
__global__ void gpu_Heat (float *h, float *g, int N) {

	// TODO: kernel computation
	//...
	
	int j = threadIdx.x + blockDim.x * blockIdx.x;
	int i = threadIdx.y + blockDim.y * blockIdx.y;
	/*if(i==0 &&  j==0){
		printf("kernel\n");
		    for (int k = 0; k <N*N; k++) {    
        printf("%f ",h[k]);
    }
	}*/

	if(i<(N-1) && j<(N-1) && i > 0 && j>0){	
		
		
		g[i*N+j]= 0.25 * (h[ i*N     + (j-1) ]+  // left
				h[ i*N     + (j+1) ]+  // right
					h[ (i-1)*N + j     ]+  // top
					h[ (i+1)*N + j     ]); // bottom
		
	
	}
}
