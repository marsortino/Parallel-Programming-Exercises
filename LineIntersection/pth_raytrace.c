#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Enabligh pthread library */
#include <pthread.h>

#include <time.h>

typedef struct{
	double x;
	double y;
	double z;
} vec;
	      

typedef struct{
	vec center;
	double r;
} block;

/* Declaring functions */
long double distance_between_line_point(vec p, vec q, vec r);
vec normal_vec(vec source, vec block_i);
long double interacting_region(vec source, vec point_on_plane, vec normal);
long double chord(double radius, double distance);
long double raytracing(int block_num, int max, block* blocklist, vec source);

/* pthread function*/
void* start_raytracing(void* size);

/* Global Variables */
int thread_count, max, **sizeslist;
long double* raytrace;
block* blocklist;
vec source;

/* pthread mutex */
pthread_mutex_t mutex;

int main(int argc, char* argv[]){

	/* Get the start time */
	struct timespec startTime, endTime;
	clock_gettime(CLOCK_MONOTONIC, &startTime);

	/* Definying variables */
	long int i = 0;
	FILE *list;
	const char *filename ="blocks_coord.txt";
	
	/* pthread objects */
	pthread_t* thread_handles;
	int thread_create_result;

	/* Get numbers of threads from command line */
	thread_count = strtol(argv[1], NULL, 10);
	thread_handles = malloc(thread_count*sizeof(pthread_t));

	/* Reading file */
	list = fopen(filename, "r");
    	if (list == NULL) {
        	fprintf(stderr, "Error: Unable to open file '%s' for reading.\n", filename);
        exit(EXIT_FAILURE);
    	}
	fscanf(list, "%d", &max);	

	/* Allocating memory*/	
	raytrace = malloc(max*sizeof(long double));
	blocklist = malloc(max*sizeof(block));
	
	/* Reading file */
	i = 0;
	while(!feof(list)){
		fscanf(list, "%lf %lf %lf %lf", &blocklist[i].center.x, &blocklist[i].center.y, &blocklist[i].center.z, &blocklist[i].r); 
		i++;
	}
	fclose(list);
	
	/* Source */ 
	source.x = 0;
	source.y = 0;
	source.z = 0;
	
	/* Defining the start/end point of each thread */	
	sizeslist = (int **) malloc(max * sizeof(int*));
	
	for (i = 0; i < max-1; i++){
		sizeslist[i] = (int*) malloc(2*sizeof(int));
		sizeslist[i][0] = i*max/thread_count;
		sizeslist[i][1] = (i+1)*max/thread_count;	
	}
	
	sizeslist[max-1] = (int*) malloc(2*sizeof(int));
	sizeslist[max-1][0] = (max-1)/thread_count;
	sizeslist[max-1][1] = max/thread_count + max % thread_count; /* In case in which max is not divided even by thread_counts */

	/* Starting the computation by each thread */
	for(i = 0; i < thread_count; i++){
		thread_create_result = pthread_create(&thread_handles[i], NULL, start_raytracing, (void*) i);
	if (thread_create_result != 0){
		fprintf(stderr, "Error: pthread_create failed with code %d\n", thread_create_result);
        	exit(EXIT_FAILURE);
   		 }
	}

	/* Freeing up memory */
	free(blocklist);
	free(sizeslist);	
	
	/* Joining the threads */
	for(i=0; i < thread_count; i++) pthread_join(thread_handles[i], NULL);

	/* Measuring end time*/
    	clock_gettime(CLOCK_MONOTONIC, &endTime);

    	/* Calculate the elapsed time */
    	double elapsedSeconds = endTime.tv_sec - startTime.tv_sec + (endTime.tv_nsec - startTime.tv_nsec) / 1e9;
	
	/* Printing the results*/
/*	for(i=0; i<max; i++){
		printf("block [%ld] raytrace value: %Le\n", i, raytrace[i]);
	}
*/	free(raytrace);

   	printf("PThreads: %d threads. Elapsed time: %f seconds\n", thread_count, elapsedSeconds);

	free(thread_handles);
	pthread_mutex_destroy(&mutex);

	return 0;
} /* main */

void* start_raytracing(void* rank){
	/* Pthread function. Each thread starts raytracing()*/
	int i, my_start, my_end, my_size;
        long my_rank = (long int) rank;

	my_start = sizeslist[my_rank][0];
	my_end = sizeslist[my_rank][1];
	my_size = my_end-my_start;
	
	long double my_raytrace[my_size];

	for(i = my_start; i < my_end; i++){
		my_raytrace[i-my_start] = raytracing(i, max, blocklist, source);
	}
	
	for(i = my_start; i < my_end; i++){
		raytrace[i] = my_raytrace[i-my_start];
	}

	return NULL;

} /* start_raytracing */
	

long double raytracing(int i, int max, block* blocklist, vec source){
	/* Defines a line passing through each sphere and the source, then it computes all the intersections */
	int j;
	long double raypath = 0, value;
	int obs_sign, block_j_sgn1, block_j_sgn2;
	vec normal;

	normal = normal_vec(source, blocklist[i].center);
	obs_sign = interacting_region(source, blocklist[i].center, normal);

	for(j=0; j < max; j++){
		if (i!=j){
			block_j_sgn1 = interacting_region(blocklist[j].center, blocklist[i].center, normal);	
			block_j_sgn2 = interacting_region(blocklist[j].center, source, normal);
			if ((obs_sign == block_j_sgn1) && (obs_sign == -block_j_sgn2)){
				value = distance_between_line_point(blocklist[i].center, source, blocklist[j].center);
				if (value == 0) raypath+=blocklist[j].r;
				else if ((value > 0) && (value < blocklist[j].r)){
					raypath += chord(blocklist[j].r, value);
				}
			}					
		}
	}
		
	return raypath;
}  /* raytracing */

vec normal_vec(vec source, vec block_i){
	/* Computes the directional vector passing through two points. */
	vec normal;

	normal.x = block_i.x-source.x;
	normal.y = block_i.y-source.y;
	normal.z = block_i.z-source.z;
	
	return normal;
} /* normal_vec */
	
long double interacting_region(vec source, vec point_on_plane, vec normal){
	/*
	 * Computes the plane passing throught the center of the i-th blob and 
	 * orthogonal to the line passing through the observer and the i-th blob center */

	long double distance;
	vec vec1;

	vec1.x = source.x-point_on_plane.x;
	vec1.y = source.y-point_on_plane.y;
	vec1.z = source.z-point_on_plane.z;
	/* Computing the dot product between vec1 and normal */
	distance = vec1.x*normal.x+vec1.y+normal.y+vec1.z+normal.z;
	distance = distance/fabsl(distance);
	return distance;
}   /* interacting_region */

long double chord(double radius, double distance){
	/* Computes the chord of the blob determined by the intersection via the blob and the rayline*/
	return 2*sqrt(2*radius*distance-distance*distance);
}  /* chord */

long double distance_between_line_point(vec p, vec q, vec r){
	/*
	 Returns the distance between a line, defined passing through two points, and a point.
        Computes the direction vector D of the line:
            D = p-q
        Computes the vector from point1 to point3.
            PR = r-p
        It then computes the projection of the vector on the direction vector:
            projection = PR x D/ D x D * D
        where 'x' stands for dot product
        Finally it normalize the distance of PR-projection
		*/
	double distance, dot1, dot2;
	vec D, PR, projection;

	D = normal_vec(p, q);

	PR.x = r.x-p.x;
	PR.y = r.y-p.y;
	PR.z = r.z-p.z;
	
	dot1 = PR.x*D.x+PR.y*D.y+PR.z*D.z;
	dot2 = D.x*D.x+D.y*D.y+D.z*D.z;
	
	projection.x = dot1/dot2*D.x;
	projection.y = dot1/dot2*D.y;
	projection.z = dot1/dot2*D.z;
	
	distance = sqrt(pow(PR.x-projection.x,2)+pow(PR.y-projection.y,2)+pow(PR.z-projection.z,2));
	return distance;
} /* distance_between_line_point */

