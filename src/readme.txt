!READ!

1. Create a job session with only 2 processes per node

srun --pty -A TG-SEE150003 -p normal -t 02:00:00 -n 4 -N 2 /bin/bash -l


This requests 4 total processes (-n) on 2 nodes (-N). Keep the # of processes twice the # of nodes for this approach.



2. Launch your application using the tacc_affinity command plus the ibrun command. And don't forget to set the # of threads per process with OMP_NUM_THREADS



3.  To start the application on all possible process:
 OMP_NUM_THREADS=8 ibrun tacc_affinity ./a.out 10mb.cypher

To start the application with determined number of process:
OMP_NUM_THREADS=16 ibrun -n 1 -o 0 tacc_affinity ./a.out 10mb.cypher  (this case is only one process)


