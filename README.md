# COMP464-Final-Project
High Performance Computing RSA Algorithm

# Compiling 

Compiling 1st version (not parallel):

<b>icpc -O3 -lgmp rsa.cpp -o rsa</b>


Compiling 2nd version (OpenMP)

<b>icpc -O3 -lgmp -fopenmp rsa.cpp -o rsa</b>


Compiling 3rd version (OpenMP + MPI)

<b>mpicc -O3 -lgmp -fopenmp rsa.cpp -o rsa</b>
