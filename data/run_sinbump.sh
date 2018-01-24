#!/usr/bin/sh

# Compute the L2-error for the sinus bump problem where the analytical solution
# is known (analytical_test=1) 
# (filling the textures with the analytical solution at each slice is
# slow compared to standard field/field distance)

# This test results can be used to verify mapping and interpolation function
# implementations

for file in *P*.json; do
    echo "../build/ffes-distance $file analytical_test=1 samples=100 | grep 'L2'"
    ../build/ffes-distance $file analytical_test=1 samples=500 | grep "L2"
    echo "-------"
done
for file in *Q*.json; do
    echo "../build/ffes-distance $file analytical_test=1 samples=100 subdiv_A=3 | grep 'L2'"
    ../build/ffes-distance $file analytical_test=1 samples=500 subdiv_A=3 | grep "L2"
    echo "-------"
done
