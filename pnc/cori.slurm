#!/bin/bash
#SBATCH -J poly2vinyl_lammps
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -N 64
#SBATCH --ntasks-per-node=64
#SBATCH -t 30:00
#SBATCH -o poly2.o%j

module load openmpi 

srun --cpu-bind=cores ~/.local/bin/lmp_custom -in polymer.lammps
