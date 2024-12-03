#include "heat.h"
#include <mpi.h>
#define min(a,b) ( ((a) < (b)) ? (a) : (b) )
#define NB 8
/*
 * Blocked Jacobi solver: one iteration step
 */
double relax_jacobi (double *u, double *utmp, unsigned sizex, unsigned sizey)
{
    double diff, sum=0.0;
    int nbx, bx, nby, by;
  
    nbx = NB;
    bx = sizex/nbx  ;
    nby = NB;
    by = sizey/nby;
    for (int ii=0; ii<nbx; ii++)
        for (int jj=0; jj<nby; jj++) 
            for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
                for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
	            utmp[i*sizey+j]= 0.25 * (u[ i*sizey     + (j-1) ]+  // left
					     u[ i*sizey     + (j+1) ]+  // right
				             u[ (i-1)*sizey + j     ]+  // top
				             u[ (i+1)*sizey + j     ]); // bottom
	            diff = utmp[i*sizey+j] - u[i*sizey + j];
	            sum += diff * diff; 
	        }

    return sum;
}

/*
 * Blocked Red-Black solver: one iteration step
 */
double relax_redblack (double *u, unsigned sizex, unsigned sizey)
{
    double unew, diff, sum=0.0;
    int nbx, bx, nby, by;
    int lsw;

    nbx = NB;
    bx = sizex/nbx;
    nby = NB;
    by = sizey/nby;
    // Computing "Red" blocks
    for (int ii=0; ii<nbx; ii++) {
        lsw = ii%2;
        for (int jj=lsw; jj<nby; jj=jj+2) 
            for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
                for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
	            unew= 0.25 * (    u[ i*sizey	+ (j-1) ]+  // left
				      u[ i*sizey	+ (j+1) ]+  // right
				      u[ (i-1)*sizey	+ j     ]+  // top
				      u[ (i+1)*sizey	+ j     ]); // bottom
	            diff = unew - u[i*sizey+ j];
	            sum += diff * diff; 
	            u[i*sizey+j]=unew;
	        }
    }

    // Computing "Black" blocks
    for (int ii=0; ii<nbx; ii++) {
        lsw = (ii+1)%2;
        for (int jj=lsw; jj<nby; jj=jj+2) 
            for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
                for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
	            unew= 0.25 * (    u[ i*sizey	+ (j-1) ]+  // left
				      u[ i*sizey	+ (j+1) ]+  // right
				      u[ (i-1)*sizey	+ j     ]+  // top
				      u[ (i+1)*sizey	+ j     ]); // bottom
	            diff = unew - u[i*sizey+ j];
	            sum += diff * diff; 
	            u[i*sizey+j]=unew;
	        }
    }

    return sum;
}

/*
 * Blocked Gauss-Seidel solver: one iteration step
 */
double relax_gauss (double *u, unsigned sizex, unsigned sizey, int rank, int numproc)
{
    printf("Una mica de guass eo no\n");
    double unew, diff, sum=0.0;
    int nbx, bx;
    MPI_Status status;

    nbx = NB;
    bx = sizex/nbx;
    //printf("\nsizex: %d, sizey %d,nbx %d, rank %d, numprocs%d \n",sizex, sizey,nbx,rank,numproc);
    //Strip mining
    //for(int i =0; i< sizex * sizey; i ++){
    //    printf(" %f",u[i]);
    //}
    
    for (int ii=0; ii<nbx; ii++){
        //printf("\nRANK %d\n",rank);
        if(rank!=0){
            MPI_Recv(&u[0], sizex, MPI_DOUBLE, rank-1,0, MPI_COMM_WORLD,&status);
            //printf("Recieved\n");
            for(int i =0; i< sizex; i ++){
                //printf(" %f",u[i]);
            }
        }

        for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) {
            for (int j=1; j<= sizey-2; j++) {
            unew= 0.25 * (    u[ i*sizey	+ (j-1) ]+  // left
                    u[ i*sizey	+ (j+1) ]+  // right
                    u[ (i-1)*sizey	+ j     ]+  // top
                    u[ (i+1)*sizey	+ j     ]); // bottom
            diff = unew - u[i*sizey+ j];
            sum += diff * diff; 
            //printf("left=%f,right= %f,top=%f,bottom= %f,res = %f ",u[ i*sizey	+ (j-1) ],u[ i*sizey	+ (j+1) ],u[ (i-1)*sizey	+ j], u[ (i+1)*sizey	+ j     ], unew);
            u[i*sizey+j]=unew;
            }
            //printf("\n");
        }
        if(rank!=numproc-1){
            //printf("Sending: \n");
            for(int i =0; i< sizex; i ++){
              //  printf(" %f",u[(sizey-2)*sizex+i]);
             }
            MPI_Send(&u[(sizey-2)*sizex], sizex, MPI_DOUBLE, rank+1,0, MPI_COMM_WORLD );
            
        }

    }
    //for(int i =0; i< sizey * bx; i ++){
    //    printf(" %f",u[i]);
    //}
    //MPI_Send(&u[0], by*sizex, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    return sum;
}

