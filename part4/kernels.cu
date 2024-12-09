#include <math.h>
#include <float.h>
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
__global__ void gpu_Heat (float *h, float *g, int N) {

	// TODO: kernel computation
	//...

	__shared__ double sh_mem[12*12];
	
	int j = threadIdx.x + blockDim.x * blockIdx.x;
	int i = threadIdx.y + blockDim.y * blockIdx.y;

	int s_idx = (threadIdx.y+1) * (blockDim.x + 2) + (threadIdx.x+1);

   

    //printf("ID:%i\n",s_idx);
    sh_mem[s_idx] = h[ i*N + j];
    if(j==0 && i==0){
		 //printf("%i\n",(blockDim.x + 2)* (blockDim.x + 2));
		for (int k = 0; k <(blockDim.x + 2)* (blockDim.x + 2); k++) {    
           // printf("%f ",sh_mem[k]);
        }
    }
    if ( threadIdx.x == 0 && j > 0) { // Left halo
        //printf("ID:%i\n",s_idx);
        sh_mem[s_idx - 1] = h[ i*N + (j-1)];
    }	
//printf("1\n");
    if ( threadIdx.x == blockDim.x - 1 && j < N - 1) { // Right halo
        sh_mem[s_idx + 1] = h[i*N  + (j+1)];
    }
    //printf("2\n");
    if ( threadIdx.y == 0 && i > 0) { // Top halo
        sh_mem[s_idx - (blockDim.x + 2)] = h[(i-1)*N + j ];
    }
       //printf("3\n");
    if ( threadIdx.y == blockDim.y - 1 && i < N - 1) { // Bottom halo
     //printf("3\n");
        sh_mem[s_idx + (blockDim.x + 2)] = h[(i+1)*N + j ];
    }
	__syncthreads();
    
    
	if(i==0 &&  j==15){
       // printf("\n huh%f\n",h[15]);
		for (int k = 0; k <(blockDim.x + 2)* (blockDim.x + 2); k++) {    
        //printf(" %f ",sh_mem[k]);
        //if(k%10==0)  printf("\n");
        }
	}
	
	


	if(i<(N-1) && j<(N-1) && i > 0 && j>0){	
		
		
		g[i*N+j]= 0.25 * (sh_mem[s_idx - 1] +  // Left
                           sh_mem[s_idx + 1] +  // Right
                           sh_mem[s_idx - (blockDim.x + 2)] +  // Top
                           sh_mem[s_idx + (blockDim.x + 2)]); // Bottom
		
	
	}
}
