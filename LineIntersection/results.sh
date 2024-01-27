#!/bin/bash

n=10
outputname="log_results.txt"
programs=("pth_raytrace" "omp_raytrace")

if [ -e "$outputname" ]; then
    rm "$outputname"
fi

touch "$outputname"

# Loop for ./raytrace
echo "Serial: raytrace" >> "$outputname"
for ((i=1; i<=n; i++))
do
    ./raytrace >> "$outputname"
done

# Loop for mpirun ./mpi_raytrace
echo "Program: mpi_raytrace" >> "$outputname"
for j in 2 4 6 8
do
    for ((i=1; i<=n; i++))
    do
        mpirun -n "$j" ./mpi_raytrace >> "$outputname"
    done
done

# Loop for other programs
for program in "${programs[@]}"
do
    echo "Program: $program" >> "$outputname"
    
    for j in 2 4 6 8
    do
        for ((i=1; i<=n; i++))
        do
            ./"$program" "$j" >> "$outputname"
        done
    done
done

bash mean.sh 

cat "results.txt"

echo "Results written to results.txt"
