#!/bin/bash
#SBATCH -J sio2_bd_20percent_lammps
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -N 16
#SBATCH -t 30:00
#SBATCH -o sio2_bd_20percent.o%j

module load openmpi 

srun --cpu-bind=cores ~/.local/bin/lmp_custom -in bd.lammps
