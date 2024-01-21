#include <stdio.h>
//#include <mpi.h>
#include <stdlib.h>

double f(double x);

int N = 100000, N_t = 100;
double a = 0, b = 3;

int main(void){
	int i;
	double* x;
	double h,y, k1, k2, k3, k4;
	x = malloc(N*sizeof(double));
	
	h = (b-a)/N;

	for(i=0; i < N; i++){
		x[i] = a+i*h;
		printf("x[%d] = %f\n", i, x[i]);
	}
	

//	int comm_sz;

//	MPI_Init();

//	MPI_Comm_size(MPI_COMM_WORLD, %comm_sz);

	for(i = 0; i < N; i++){
		y = f(x[i]);
		k1 = f(y);
		k2 = f(y+h*k1/2);
		k3 = f(y+h*k2/2);
		k4 = f(y+h*k3);
		printf("viene: %f\n",h/6*(k1+2*k2+2*k3+k4));
		x[i] = x[i] + h/6*(k1+2*k2+2*k3+k4);
		printf("x[%d] = %f\n", i, x[i]);
	}		

//	MPI_Finalize();
	return 0;
}

double f(double x){
	/* Compute the function 3x*2 */
	return 3*x*x;
}


