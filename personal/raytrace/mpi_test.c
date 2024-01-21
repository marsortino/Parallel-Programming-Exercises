#include <stdio.h>
#include <mpi.h>

typedef struct {
    int x;
    int y;
} Point;

typedef struct {
    int id;
    Point coordinates;
} Particle;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define MPI datatypes for nested structs
    MPI_Datatype pointType, particleType;
    MPI_Aint offsets[2];
    int blocklengths[2];
    MPI_Datatype types[2];

    // Define the Point datatype
    offsets[0] = 0;
    types[0] = MPI_INT;
    blocklengths[0] = 2;
    MPI_Type_create_struct(1, blocklengths, offsets, types, &pointType);
    MPI_Type_commit(&pointType);

    // Define the Particle datatype
    offsets[0] = 0;
    offsets[1] = sizeof(int);
    types[1] = pointType;
    blocklengths[1] = 1;
    MPI_Type_create_struct(2, blocklengths, offsets, types, &particleType);
    MPI_Type_commit(&particleType);

    // Example data
    Particle particle;
    particle.id = rank;
    particle.coordinates.x = 10 * rank;
    particle.coordinates.y = 20 * rank;

    // Broadcast the Particle struct using MPI datatype
    MPI_Bcast(&particle, 1, particleType, 0, MPI_COMM_WORLD);

    // Print received data
    printf("Rank %d received Particle: ID=%d, Coordinates=(%d, %d)\n",
           rank, particle.id, particle.coordinates.x, particle.coordinates.y);

    MPI_Type_free(&pointType);
    MPI_Type_free(&particleType);
    MPI_Finalize();
    return 0;
}

