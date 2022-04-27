#!/bin/bash
#SBATCH -J test_lammps
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -t 30:00
#SBATCH -o test_lammps.o%j

module load lammps

srun -n 4 -c 1 --cpu-bind=cores lmp_cori -in input.lammps
