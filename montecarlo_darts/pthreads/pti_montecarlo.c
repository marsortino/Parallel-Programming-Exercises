#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>

int thread_count;
double number_in_circle = 0;
pthread_mutex_t mutex;
void* pi_estimate(void* number_of_tosses);

int main(int argc, char* argv[]){
	int thread_count;
	long thread;
	pthread_t* thread_handles;
	long long int number_of_tosses, number_of_tosses_per_thread, remainder;

	thread_count = strtol(argv[1], NULL, 10);

	/* Get numbers of threads from command line */
	thread_handles = malloc(thread_count*sizeof(pthread_t));

	printf("How many darts do you want to throw?\n");
	scanf("%lld", &number_of_tosses);

	number_of_tosses_per_thread = number_of_tosses / thread_count;
	remainder = number_of_tosses_per_thread % thread_count;
	printf("number of tosses per thread %lld\n", number_of_tosses_per_thread);

	pthread_mutex_init(&mutex, NULL);

	/* Starting the first thread with n_tosses + remainder */	
	pthread_create(&thread_handles[0], NULL, pi_estimate, (void*) number_of_tosses_per_thread + remainder);
	/* Starting the others threads */
	for (thread = 1; thread < thread_count; thread++){
		pthread_create(&thread_handles[thread], NULL, pi_estimate, (void*) number_of_tosses_per_thread);
	}
	/* Joining the threads */
	for (thread = 0; thread < thread_count; thread++){
		pthread_join(thread_handles[thread], NULL);
	}
	
	free(thread_handles);
	pthread_mutex_destroy(&mutex);

	/* Printing the final result */
	printf("pi = %f\n", 4*number_in_circle/number_of_tosses);
	return 0;
}

void* pi_estimate(void* n_tosses){
	double x,y, distance_squared;
	long long int i;

 	long long int my_n_tosses = (long long int) n_tosses;
	printf("starting\n");
	for (i = 0; i < my_n_tosses; i++){
		x = rand() % RAND_MAX;
		x = x / RAND_MAX*2 - 1;
		y = rand() % RAND_MAX;
		y = y / RAND_MAX*2 - 1;

		distance_squared = x*x+y*y;
		if (distance_squared <= 1){
			pthread_mutex_lock(&mutex);       
			number_in_circle++;
			pthread_mutex_unlock(&mutex);
		}
	}
	return NULL;
}

