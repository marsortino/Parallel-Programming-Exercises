#!/bin/bash

# Check if results.txt exists, and delete it if it does
if [ -e results.txt ]; then
    rm results.txt
fi

# Extract elapsed times from the provided data in log_results.txt
serial_times=$(grep "Serial version" log_results.txt | awk '{print $5}')
mpi_2_times=$(grep "MPI: 2 processors" log_results.txt | awk '{print $6}')
mpi_4_times=$(grep "MPI: 4 processors" log_results.txt | awk '{print $6}')
mpi_6_times=$(grep "MPI: 6 processors" log_results.txt | awk '{print $6}')
mpi_8_times=$(grep "MPI: 8 processors" log_results.txt | awk '{print $6}')
pthread_2_times=$(grep "PThreads: 2 threads" log_results.txt | awk '{print $6}')
pthread_4_times=$(grep "PThreads: 4 threads" log_results.txt | awk '{print $6}')
pthread_6_times=$(grep "PThreads: 6 threads" log_results.txt | awk '{print $6}')
pthread_8_times=$(grep "PThreads: 8 threads" log_results.txt | awk '{print $6}')
omp_2_times=$(grep "OMP: 2 threads" log_results.txt | awk '{print $6}')
omp_4_times=$(grep "OMP: 4 threads" log_results.txt | awk '{print $6}')
omp_6_times=$(grep "OMP: 6 threads" log_results.txt | awk '{print $6}')
omp_8_times=$(grep "OMP: 8 threads" log_results.txt | awk '{print $6}')

# Function to calculate average
calculate_average() {
  local sum=0
  local count=0
  for time in $1; do
    sum=$(echo "$sum + $time" | bc)
    count=$((count + 1))
  done
  echo "scale=3; ($sum) / $count" | bc
}

# Function to format elapsed time
format_elapsed_time() {
  local time=$1
  if [ $(echo "$time < 1" | bc) -eq 1 ]; then
    echo "0$time"
  else
    echo "$time"
  fi
}

# Calculate the average for each subset
average_serial=$(format_elapsed_time "$(calculate_average "$serial_times")")
average_mpi_2=$(format_elapsed_time "$(calculate_average "$mpi_2_times")")
average_mpi_4=$(format_elapsed_time "$(calculate_average "$mpi_4_times")")
average_mpi_6=$(format_elapsed_time "$(calculate_average "$mpi_6_times")")
average_mpi_8=$(format_elapsed_time "$(calculate_average "$mpi_8_times")")
average_pthread_2=$(format_elapsed_time "$(calculate_average "$pthread_2_times")")
average_pthread_4=$(format_elapsed_time "$(calculate_average "$pthread_4_times")")
average_pthread_6=$(format_elapsed_time "$(calculate_average "$pthread_6_times")")
average_pthread_8=$(format_elapsed_time "$(calculate_average "$pthread_8_times")")
average_omp_2=$(format_elapsed_time "$(calculate_average "$omp_2_times")")
average_omp_4=$(format_elapsed_time "$(calculate_average "$omp_4_times")")
average_omp_6=$(format_elapsed_time "$(calculate_average "$omp_6_times")")
average_omp_8=$(format_elapsed_time "$(calculate_average "$omp_8_times")")

# Write the results to results.txt
echo "Average elapsed time for each subset:" >> results.txt
echo "Serial version: $average_serial seconds" >> results.txt
echo "MPI: 2 processors - $average_mpi_2 seconds" >> results.txt
echo "MPI: 4 processors - $average_mpi_4 seconds" >> results.txt
echo "MPI: 6 processors - $average_mpi_6 seconds" >> results.txt
echo "MPI: 8 processors - $average_mpi_8 seconds" >> results.txt
echo "PThreads: 2 threads - $average_pthread_2 seconds" >> results.txt
echo "PThreads: 4 threads - $average_pthread_4 seconds" >> results.txt
echo "PThreads: 6 threads - $average_pthread_6 seconds" >> results.txt
echo "PThreads: 8 threads - $average_pthread_8 seconds" >> results.txt
echo "OMP: 2 threads - $average_omp_2 seconds" >> results.txt
echo "OMP: 4 threads - $average_omp_4 seconds" >> results.txt
echo "OMP: 6 threads - $average_omp_6 seconds" >> results.txt
echo "OMP: 8 threads - $average_omp_8 seconds" >> results.txt

