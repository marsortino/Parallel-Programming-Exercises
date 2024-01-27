#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

double montecarlo(long long int);

void main(void){
	int comm_sz;
	int my_rank;
	
	long long int local_n;
	double local_result, global_result, pi_estimate;
	MPI_Init(NULL, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (my_rank == 0){
		long long int remainder, number_of_tosses = 100000000;
		
		local_n = number_of_tosses/comm_sz;
		remainder = number_of_tosses % comm_sz;
		MPI_Bcast(&local_n, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		local_result = montecarlo(local_n + remainder);
		MPI_Reduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);	
		pi_estimate = 4*global_result/((double) number_of_tosses);
		printf("result: %f\n", pi_estimate);
	}
	else{
		MPI_Bcast(&local_n, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
		local_result = montecarlo(local_n);
		MPI_Reduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);	
	}
	MPI_Finalize();
}

double montecarlo(long long int number_of_tosses){
	long long int number_in_circle, toss;
	double x,y, distance_squared, pi_estimate;
	
	number_in_circle = 0;

	for (toss = 0; toss < number_of_tosses; toss++){
	
		x = rand() % RAND_MAX;
		y = rand() % RAND_MAX;
		x = x/RAND_MAX*2-1;
		y = y/RAND_MAX*2-1;
	
		distance_squared = x*x + y*y;
		if (distance_squared <= 1) number_in_circle++;
	}
	
	return number_in_circle;
}
