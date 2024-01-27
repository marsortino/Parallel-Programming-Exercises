#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#	include <omp.h>
#endif

int main(int argc, char* argv[]){

	int i, thread_count = strtol(argv[1], NULL, 10), number_of_toss, number_in_circle = 0;
	double x, y, estimate, distance;
	printf("number of tosses?\n");
	scanf("%d", &number_of_toss);
# 	pragma omp parallel for num_threads(thread_count) reduction(+:number_in_circle) private(x,y,distance)
	for (i=0; i < number_of_toss; i++){
		
		x = rand() % RAND_MAX;
		y = rand() % RAND_MAX;
		x = x / RAND_MAX * 2 - 1;
		y = y / RAND_MAX * 2 - 1;

		distance = x*x+y*y;
		if (distance <= 1.0) number_in_circle++;
	}
	estimate = ((double) number_in_circle / number_of_toss) * 4;
	printf("pi is: %f\n", estimate);
	return 0;
}
