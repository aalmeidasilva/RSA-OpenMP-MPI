# COMP464-Final-Project
High Performance Computing RSA Algorithm

# Compiling Instructions

Compiling 1st version (not parallel)
# icpc -O3 -lgmp rsa.cpp -o rsa

Compiling 2nd version (OpenMP)
# icpc -O3 -lgmp -fopenmp rsa.cpp -o rsa

Compiling 3rd version (OpenMP + MPI)
# mpicc -O3 -lgmp -fopenmp rsa.cpp -o rsa
