#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>

/* Enabling MPI */
#include <mpi.h>

typedef struct{
	double x;
	double y;
	double z;
} vec;
	      

typedef struct{
	vec center;
	double r;
} block;

/* Functions declaration */ 
long double distance_between_line_point(vec p, vec q, vec r);
vec normal_vec(vec source, vec block_i);
long double interacting_region(vec source, vec point_on_plane, vec normal);
long double chord(double radius, double distance);
long double raytracing(int block_num, int max, block* blocklist, vec source);

/* MPI functions declaration */
void Build_mpi_type_block(MPI_Datatype* input_block, MPI_Datatype input_vec);
void Build_mpi_type_vec(MPI_Datatype* input_vec);
/* */

/* MAIN */
int main(void){

	MPI_Init(NULL, NULL);
	int my_rank, comm_sz, local_max, local_i;
	
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	double local_start, local_finish, local_elapsed, elapsed;
	vec source;
	int i = 0, max = 2824, tmp_local_max, remainder;
	FILE *list;
	block *blocklist;

	
	int local_max_counts[comm_sz], local_displs[comm_sz];
	
	MPI_Datatype MPI_vec, MPI_block;
	Build_mpi_type_vec(&MPI_vec);
	Build_mpi_type_block(&MPI_block, MPI_vec);

	if (my_rank == 0){
		
		list = fopen("blocks_coord.txt", "r");
/*    		if (list == NULL) {
      		  	fprintf(stderr, "Error opening file\n");
       			 MPI_Abort(MPI_COMM_WORLD, 1);
   		 }
*/
		fscanf(list, "%d", &max);
		blocklist = malloc((32*max)*sizeof(MPI_INT));
	       	for(i=0; i < max; i++){
        	        fscanf(list, "%lf %lf %lf %lf", &blocklist[i].center.x, &blocklist[i].center.y, &blocklist[i].center.z, &blocklist[i].r);
	       	}
		fclose(list);
	}
	
	/* Starting measuring time */
	MPI_Barrier(MPI_COMM_WORLD);
	local_start = MPI_Wtime();
	
	/* Starting Bcasting data to each processor */
	MPI_Bcast(&max, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (my_rank != 0) blocklist = malloc((32*max)*sizeof(MPI_INT));
	
	MPI_Bcast(blocklist, max, MPI_block, 0, MPI_COMM_WORLD);
	
	/* Starting raytracing procedures */
	tmp_local_max = max/comm_sz;
	local_max = (my_rank+1)*max/comm_sz;
	local_i = my_rank*max/comm_sz;
	remainder = max % comm_sz;
	
	if (my_rank == comm_sz-1){
			local_max += remainder;
			tmp_local_max += remainder;
		}

	/* Raytracing */	
	long double *global_raytrace, local_raytrace[tmp_local_max];
	
	source.x = 0;
	source.y = 0;
	source.z = 0;
	for ( i = local_i; i<local_max; i++){
		local_raytrace[i-local_i] = raytracing(i, max, blocklist, source);
	}
	
	free(blocklist);
	/* MPI_Gatherv procedures */ 
	global_raytrace = malloc(max*sizeof(long double));
		
	local_displs[0] = 0;
	if (my_rank == comm_sz-1) tmp_local_max = max/comm_sz;	
	for (i=0; i<comm_sz; i++){
		local_max_counts[i] = tmp_local_max;
		if (i<comm_sz-1) local_displs[i+1] = local_displs[i]+local_max_counts[i];
	}

	local_max_counts[comm_sz-1] = tmp_local_max+remainder;

	MPI_Gatherv(local_raytrace, local_max_counts[my_rank], MPI_LONG_DOUBLE, global_raytrace, local_max_counts, local_displs, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
	
	/* Measuring end time */
	MPI_Barrier(MPI_COMM_WORLD);
	local_finish = MPI_Wtime();
	local_elapsed = local_finish-local_start;
	MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (my_rank == 0){
		 //for (i=0; i<max; i++) printf("raytracing[%d]: %Le\n", i, global_raytrace[i]);
		printf("MPI: %d processors. Elapsed time: %f seconds\n", comm_sz, elapsed);
	}
	free(global_raytrace);
	MPI_Type_free(&MPI_block);
	MPI_Type_free(&MPI_vec);
	MPI_Finalize();
	return 0;
}
	
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
}  /* normal_vec */	
	
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
	long double distance, dot1, dot2;
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
}  /* distance_between_line_point */

void Build_mpi_type_vec(MPI_Datatype* input_vec){
	/* Define a MPI Datatype with structure vec (i.e. double x double y double z) */
	int array_of_blocklenghts_vec[3] = {1, 1, 1};
	MPI_Datatype array_of_types_vec[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Aint disps_vec[3] = {offsetof(vec,x), offsetof(vec,y), offsetof(vec,z)};

	MPI_Type_create_struct(3, array_of_blocklenghts_vec, disps_vec, array_of_types_vec, input_vec);
	MPI_Type_commit(input_vec);
}   /* Build_mpi_type_vec */ 

void Build_mpi_type_block(MPI_Datatype* input_block, MPI_Datatype input_vec){
	/* Define a MPI Datatype with block structure (i.e. vec vector, double radius) */
	int array_of_blocklenghts_block[2] = {1, 1};
	MPI_Datatype array_of_types_block[2] = {input_vec, MPI_DOUBLE};
 	MPI_Aint disps_block[2] = {offsetof(block, center), offsetof(block, r)};
	
	MPI_Type_create_struct(2, array_of_blocklenghts_block, disps_block, array_of_types_block, input_block);
	MPI_Type_commit(input_block);
}   /* Build_mpi_type_block */
