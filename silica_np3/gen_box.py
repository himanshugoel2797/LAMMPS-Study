#Generate file containing 250nm atoms in grid in lammps format
import os
import numpy as np
import sys

out_grid_side = 20e-6
grid_side = out_grid_side - 250.e-09
atom_list = []
atom_count = 0
tgt_atom_count = 156455
atoms_per_side = int(tgt_atom_count ** (1./3) + 1)
atom_spacing = (2 * grid_side) / atoms_per_side

def tform(x):
    return x * atom_spacing - grid_side

def gen_grid():
    global atom_count
    global atom_list

    for i in range(atoms_per_side):
        for j in range(atoms_per_side):
            for k in range(atoms_per_side):
                atom_list.append([tform(i),tform(j),tform(k)])
                atom_count += 1
                if atom_count == tgt_atom_count:
                    return

gen_grid()

#Write file
with open('grid.lammps','w') as f:
    f.write('LAMMPS data file for grid\n')
    f.write('\n')
    f.write('{} atoms\n'.format(atom_count))
    f.write('{} atom types\n'.format(1))
    f.write('\n')
    f.write('{} {} xlo xhi\n'.format(-out_grid_side,out_grid_side))
    f.write('{} {} ylo yhi\n'.format(-out_grid_side,out_grid_side))
    f.write('{} {} zlo zhi\n'.format(-out_grid_side,out_grid_side))
    f.write('\n')
    f.write('Atoms\n')
    f.write('\n')
    for i in range(atom_count):
        f.write('{} {} 125.e-9 2650 {} {} {}\n'.format(i+1,1,atom_list[i][0],atom_list[i][1],atom_list[i][2]))
    