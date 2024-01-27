#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

long double distance_between_line_point(vec p, vec q, vec r);
vec normal_vec(vec source, vec block_i);
long double interacting_region(vec source, vec point_on_plane, vec normal);
long double chord(double radius, double distance);
long double raytracing(int block_num, int max, block* blocklist, vec source);

int main(void){
	struct timespec start, end;
	double elapsed_time;
	
	
	clock_gettime(CLOCK_MONOTONIC, &start);


	vec source;
	int i = 0, max;
	FILE *list;
	block *blocklist;
	long double *raytrace;
	
	source.x = 0;
	source.y = 0;
	source.z = 0;

	list = fopen("blocks_coord.txt", "r");
	fscanf(list, "%d", &max);	
	blocklist = malloc(max*sizeof(block));
	while(!feof(list)){
		fscanf(list, "%lf %lf %lf %lf", &blocklist[i].center.x, &blocklist[i].center.y, &blocklist[i].center.z, &blocklist[i].r); 
		i++;
	}
	fclose(list);
	
	raytrace = malloc(max*sizeof(block));

	for(i=0; i<max; i++){
		raytrace[i] = raytracing(i, max, blocklist, source);
	//	printf("block [%d] raytrace value: %Le\n", i, raytrace[i]);
	}
	
	clock_gettime(CLOCK_MONOTONIC, &end);
	elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        printf("Serial version. Elapsed time: %f seconds\n", elapsed_time);
	return 0;
}
	
long double raytracing(int i, int max, block* blocklist, vec source){
	/* Fare il computo di raytracing come nel formato python. block_num è il numero del blocco da computare. max è la dimensione dell'array di x y e z. */
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
}	
	
long double interacting_region(vec source, vec point_on_plane, vec normal){
	/*
	 * Computes the plane passing throught the center of the i-th blob and 
	 * orthogonal to the line passing through the observer and the i-th blob center */

	double distance;
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
}



