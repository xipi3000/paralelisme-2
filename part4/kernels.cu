#include <math.h>
#include <float.h>
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
__global__ void gpu_Heat (float *h, float *g, int N) {

	// TODO: kernel computation
	//...

	extern  __shared__ double sh_mem[];
	
	int j = threadIdx.x + blockDim.x * blockIdx.x;
	int i = threadIdx.y + blockDim.y * blockIdx.y;

	int s_idx = threadIdx.y * (blockDim.x + 2) + threadIdx.x;


	
    if ( threadIdx.x == 0 && j > 0) { // Left halo
        sh_mem[s_idx - 1] = h[ i*N + (j-1)];
    }	
	printf("kernel\n");
    if ( threadIdx.x == blockDim.x - 1 && j < N - 1) { // Right halo
        sh_mem[s_idx + 1] = h[i*N  + (j+1)];
    }
    if ( threadIdx.y == 0 && i > 0) { // Top halo
        sh_mem[s_idx - (blockDim.x + 2)] = h[(i-1)*N + j ];
    }
    if ( threadIdx.y == blockDim.y - 1 && i < N - 1) { // Bottom halo
        sh_mem[s_idx + (blockDim.x + 2)] = h[(i+1)*N + j ];
    }
	__syncthreads();

	if(i==0 &&  j==0){
		printf("kernel\n");
		    for (int k = 0; k <N*N; k++) {    
        printf("%f ",sh_mem[k]);
    }
	}
	
	


	if(i<(N-1) && j<(N-1) && i > 0 && j>0){	
		
		
		g[i*N+j]= 0.25 * (sh_mem[s_idx - 1] +  // Left
                           sh_mem[s_idx + 1] +  // Right
                           sh_mem[s_idx - (blockDim.x + 2)] +  // Top
                           sh_mem[s_idx + (blockDim.x + 2)]); // Bottom
		
	
	}
}
