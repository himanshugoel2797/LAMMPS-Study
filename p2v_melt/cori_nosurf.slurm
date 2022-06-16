#!/bin/bash
#SBATCH -J sio2_pnc_lammps
#SBATCH -C haswell
#SBATCH -q debug  
#SBATCH -N 32
#SBATCH --ntasks-per-node=32
#SBATCH -t 30:00
#SBATCH -o sio2_pnc_%j

module load openmpi 

srun --cpu-bind=cores ~/.local/bin/lmp_custom2 -in bd_nosurf.lammps -sf opt 
