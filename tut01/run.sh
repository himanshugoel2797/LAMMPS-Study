#!/bin/sh
lmp_serial -in input.lammps
#mpirun -n 4 lmp_mpi -in input.lammps
