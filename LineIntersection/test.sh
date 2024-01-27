#!/bin/bash

n=5  # Number of iterations
total=0

for ((i=1; i<=n; i++))
do
    # Replace this line with the actual command that generates the output
    output=$(./raytrace | awk '/Elapsed time:/ {print $(NF-1)}')

    # Convert scientific notation to a float
    float_value=$(printf "%.6f" "$output")

    # Add the float value to the total
    total=$(awk "BEGIN {printf \"%.6f\", $total + $float_value}")

    # Print the current float value (optional)
    echo "Iteration $i: $float_value"
done

# Calculate the mean
mean=$(awk "BEGIN {printf \"%.6f\", $total / $n}")

# Print the mean
echo "Mean Elapsed Time: $mean"

