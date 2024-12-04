/*
 * Iterative solver for heat distribution
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "heat.h"

void usage( char *s )
{
    fprintf(stderr, 
	    "Usage: %s <input file> [result file]\n\n", s);
}

int main( int argc, char *argv[] )
{
    unsigned iter;
    FILE *infile, *resfile;
    char *resfilename;
    int myid, numprocs;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

if (myid == 0) {
      printf("I am the master (%d) and going to distribute work to %d additional workers ...\n", myid, numprocs-1);

    // algorithmic parameters
    algoparam_t param;
    int np;

    double runtime, flop;
    double residual=0.0;

    // check arguments
    if( argc < 2 )
    {
	usage( argv[0] );
	return 1;
    }

    // check input file
    if( !(infile=fopen(argv[1], "r"))  ) 
    {
	fprintf(stderr, 
		"\nError: Cannot open \"%s\" for reading.\n\n", argv[1]);
      
	usage(argv[0]);
	return 1;
    }

    // check result file
    resfilename= (argc>=3) ? argv[2]:"heat.ppm";

    if( !(resfile=fopen(resfilename, "w")) )
    {
	fprintf(stderr, 
		"\nError: Cannot open \"%s\" for writing.\n\n", 
		resfilename);
	usage(argv[0]);
	return 1;
    }

    // check input
    if( !read_input(infile, &param) )
    {
	fprintf(stderr, "\nError: Error parsing input file.\n\n");
	usage(argv[0]);
	return 1;
    }
    print_params(&param);

    // set the visualization resolution
    
    param.u     = 0;
    param.uhelp = 0;
    param.uvis  = 0;
    param.visres = param.resolution;
   
    if( !initialize(&param) )
	{
	    fprintf(stderr, "Error in Solver initialization.\n\n");
	    usage(argv[0]);
            return 1;
	}


    // full size (param.resolution are only the inner points)
    np = param.resolution + 2;
    unsigned num_rows = (param.resolution / numprocs) +2;
    unsigned num_rows_to_compute = (param.resolution / numprocs);
    // starting time
    runtime = wtime();
    unsigned last_rows = (param.resolution - num_rows_to_compute * (numprocs-1))+2;
    //printf("num_rows= %u, rows_comp = %u ,last_rows = %u,numprocs = %i\n",num_rows,num_rows_to_compute,last_rows,numprocs);
    //for (int i = 0; i <(sizeof(param.u)*sizeof(param.u)); i++) {
    //param.u[i] = (double)(i + 1);
    //}

    // send to workers the necessary data to perform computation
    for (int i=1; i<numprocs; i++) {
        MPI_Send(&param.maxiter, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        MPI_Send(&np, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        if(i==numprocs-1){
            //printf("Enviant last rows\n");
            MPI_Send(&last_rows, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Send(&num_rows, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        MPI_Send(&param.algorithm, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        //printf("Master- Num columns = %d\n",np );
        //printf("Master- Num rows = %d\n",num_rows );
        if(i==numprocs-1){
            //printf("Enviant matrix desde %u fins %u\n",(i*(num_rows_to_compute* np)),((i*(num_rows_to_compute* np))+(last_rows * np)));
            MPI_Send(&param.u[i*(num_rows_to_compute* np)],(last_rows * np), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }else
        {
        MPI_Send(&param.u[i*(num_rows_to_compute* np)], (num_rows * np), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        //MPI_Send(&param.uhelp[0], (np)*(np), MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
        //printf("Sending work to workers...\n");
       
	
    }


    iter = 0;
    while(1) {
	switch( param.algorithm ) {
	    case 0: // JACOBI
	            residual = relax_jacobi(param.u, param.uhelp, np, np);
		    // Copy uhelp into u
		    for (int i=0; i<np; i++)
    		        for (int j=0; j<np; j++)
	    		    param.u[ i*np+j ] = param.uhelp[ i*np+j ];
		    break;
	    case 1: // RED-BLACK
		    residual = relax_redblack(param.u, np, np);
		    break;
	    case 2: // GAUSS
		    residual = relax_gauss(param.u, np, num_rows,myid,numprocs);
		    break;
	    }

        iter++;

        // solution good enough ?
        if (residual < 0.00005) break;

        // max. iteration reached ? (no limit with maxiter=0)
        if (param.maxiter>0 && iter>=param.maxiter) break;
    }
    //for (int i = (num_rows_to_compute+1)*np; i <(sizeof(param.u)*sizeof(param.u)); i++) {
    //    param.u[i] = (double)(1);
     //   printf("%f ",param.u[i]);
    //}
    int next_cell= (num_rows_to_compute+1)* np;
    for (int i=1; i<numprocs-1; i++) {
        //printf("\nRecv from: %i\n ",i);
        MPI_Recv(&param.u[next_cell], (num_rows_to_compute* np), MPI_DOUBLE, i,0, MPI_COMM_WORLD,&status );
        // printf("\nNext cell : %i\n ",next_cell);
        next_cell += num_rows_to_compute* np;
    }
    MPI_Recv(&param.u[next_cell], ((last_rows-2)* np), MPI_DOUBLE, numprocs-1,0, MPI_COMM_WORLD,&status );


   // MPI_Reduce(&residual, &residual, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

    //for (int i = 0; i <(sizeof(param.u)*sizeof(param.u)); i++) {
        
    //    printf("%f ",param.u[i]);
    //}
    // Flop count after iter iterations
    flop = iter * 11.0 * param.resolution * param.resolution;
    // stopping time
    runtime = wtime() - runtime;

    fprintf(stdout, "Time: %04.3f ", runtime);
    fprintf(stdout, "(%3.3f GFlop => %6.2f MFlop/s)\n", 
	    flop/1000000000.0,
	    flop/runtime/1000000);
    fprintf(stdout, "Convergence to residual=%f: %d iterations\n", residual, iter);

    // for plot...
    coarsen( param.u, np, np,
	     param.uvis, param.visres+2, param.visres+2 );
  
    write_image( resfile, param.uvis,  
		 param.visres+2, 
		 param.visres+2 );

    finalize( &param );

    fprintf(stdout, "Process %d finished computing with residual value = %f\n", myid, residual);

    MPI_Finalize();

    return 0;

} else {

    printf("I am worker %d and ready to receive work to do ...\n", myid);

    // receive information from master to perform computation locally

    int columns, rows, np;
    int iter, maxiter;
    int algorithm;
    double residual =0;

    MPI_Recv(&maxiter, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&columns, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&algorithm, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    
    //printf("Worker %i => Rows = %d, Columns = %d\n",myid,rows, columns);
    // allocate memory for worker
    double * u = calloc( sizeof(double),((rows*columns)));
    double * uhelp = calloc( sizeof(double),(rows+2)*(columns+2) );

    if( (!u) || (!uhelp) )
    {
        fprintf(stderr, "Error: Cannot allocate memory\n");
        return 0;
    }
    
    // fill initial values for matrix with values received from master
    //printf("RANG = %i Waiting  ...\n",myid);
    MPI_Recv(&u[0], (rows*columns), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    //for(int i =0; i< rows * columns; i ++){
    //    printf(" %f",u[i]);
    //}
    //MPI_Recv(&uhelp[0], (rows+2)*(columns+2), MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status);
    //printf("Work recieved again -rank %i...\n",myid);
    iter = 0;
    while(1) {
	switch( algorithm ) {
	    case 0: // JACOBI
	        residual = relax_jacobi(u, uhelp, np, np);
		    // Copy uhelp into u
		    for (int i=0; i<np; i++)
    		        for (int j=0; j<np; j++)
	    		    u[ i*np+j ] = uhelp[ i*np+j ];
		    break;
	    case 1: // RED-BLACK
		    residual = relax_redblack(u, np, np);
		    break;
	    case 2: // GAUSS
		    residual = relax_gauss(u, columns, rows,myid,numprocs);
		    break;
	    }

        iter++;

        // solution good enough ?
        if (residual < 0.00005){
            printf("Residual\n");
             break;}

        // max. iteration reached ? (no limit with maxiter=0)
        if (maxiter>0 && iter>=maxiter) {
            printf("maxiter\n");
            break;}
    }
    MPI_Send(&u[columns], (rows-2)*columns, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    //double global_sum;
    //MPI_Reduce(&residual, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
    if( u ) free(u); 
    if( uhelp ) free(uhelp);

    fprintf(stdout, "Process %d finished computing %d iterations with residual value = %f\n", myid, iter, residual);

    MPI_Finalize();
    exit(0);
  }
}
