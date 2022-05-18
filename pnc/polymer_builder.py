# Copyright (c) 2022 Himanshu Goel
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

import os
import sys
import numpy as np

chain_length = 86 # corresponds to 9kg/mol
one_d = 1#5.29 # Angstrom
H_d = 0.3 #d units
H_dist = 0.37 #d units
H_mass = 17.007 #g/mol for OH group
B_mass = 105.14 - H_mass #g/mol for the r est of the monomer
box_side = 200 #529nm^3

atoms = []
bonds = []
angles = []
dihedrals = []

max_x = 0
max_y = 0
max_z = 0

def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(npoints, ndim)
    vec /= np.linalg.norm(vec, axis=1, keepdims=True)
    return vec

def perpendicular_vector(v):
    if v[1] == 0 and v[2] == 0:
        if v[0] == 0:
            raise ValueError('zero vector')
        else:
            return np.cross(v, [0, 1, 0])
    return np.cross(v, [1, 0, 0])

def build_monomer(chain_idx, monomer_idx, origin, direction):
    global max_x
    global max_y
    global max_z
    b_id = len(atoms)
    h_id = len(atoms) + 1

    #generate random 3d vector using numpy
    b_pos = origin + direction * monomer_idx
    h_pos = b_pos + perpendicular_vector(direction) * H_dist

    atoms.append((b_id, chain_idx, 1, b_pos))
    atoms.append((h_id, chain_idx, 2, h_pos)) #h bead is 90 degrees from b bead

    if max_x < b_pos[0]:
        max_x = b_pos[0]
    if max_y < b_pos[1]:
        max_y = b_pos[1]
    if max_z < b_pos[2]:
        max_z = b_pos[2]

    if max_x < h_pos[0]:
        max_x = h_pos[0]
    if max_y < h_pos[1]:
        max_y = h_pos[1]
    if max_z < h_pos[2]:
        max_z = h_pos[2]
        
    #Bonds
    if monomer_idx > 0: #bond to previous monomer
        bonds.append((len(bonds), 1, b_id - 2, b_id))
    bonds.append((len(bonds), 2, b_id, h_id)) #bond to h bead

    #Angles
    if monomer_idx > 0: #angle to previous monomer
        angles.append((len(angles), 1, b_id - 2, b_id, h_id))

def build_chain(chain_idx, chain_len, origin, direction):
    for i in range(chain_len):
        build_monomer(chain_idx, i, origin, direction)

def add_nanoparticle(pos, rad, chain_cnt, chain_len, chain_idx):
    np_idx = len(atoms)
    atoms.append((np_idx, chain_idx, 3, pos))
    
    origins = sample_spherical(chain_cnt, 3)
    for i in range(chain_cnt):
        bonds.append((len(bonds), 3, np_idx, len(atoms)))
        build_chain(chain_idx, chain_len, pos + origins[i] * rad, origins[i])


def output_polymer(filename='polymer.txt'):
    #open output file
    out_file = open(filename, 'w')
    out_file.write('# Model for poly(2-vinyl) chain\n')
    out_file.write('\n')
    out_file.write('{} atoms\n'.format(len(atoms)))
    out_file.write('{} bonds\n'.format(len(bonds)))
    out_file.write('{} angles\n'.format(len(angles)))
    out_file.write('\n')
    out_file.write('\n')
    out_file.write('{} atom types\n'.format(3))
    out_file.write('{} bond types\n'.format(3))
    out_file.write('{} angle types\n'.format(1))
    out_file.write('\n')
    out_file.write('{} {} xlo xhi\n'.format(-200, 200))
    out_file.write('{} {} ylo yhi\n'.format(-200, 200))
    out_file.write('{} {} zlo zhi\n'.format(-200, 200))
    out_file.write('\n')
    out_file.write('\n')
    out_file.write('Masses\n')
    out_file.write('\n')
    out_file.write('1 {}\n'.format(B_mass))
    out_file.write('2 {}\n'.format(H_mass))
    out_file.write('3 {}\n'.format(60.06))
    out_file.write('\n')
    out_file.write('\n')
    out_file.write('Atoms\n')
    out_file.write('\n')
    box_half = box_side / 2
    for atom_id, chain_idx, atom_type, pos in atoms:
        out_file.write('{} {} {} {} {} {}\n'.format(atom_id + 1, chain_idx + 1, atom_type, (pos[0] + 0.5) * one_d, (pos[1] + 0.5) * one_d, (pos[2] + 0.5) * one_d))
    out_file.write('\n')
    out_file.write('\n')
    out_file.write('Bonds\n')
    out_file.write('\n')
    for bond_id, bond_type, atom_id_1, atom_id_2 in bonds:
        out_file.write('{} {} {} {}\n'.format(bond_id + 1, bond_type, atom_id_1 + 1, atom_id_2 + 1))
    out_file.write('\n')
    out_file.write('\n')
    out_file.write('Angles\n')
    out_file.write('\n')
    for angle_id, angle_type, atom_id_1, atom_id_2, atom_id_3 in angles:
        out_file.write('{} {} {} {} {}\n'.format(angle_id + 1, angle_type, atom_id_1 + 1, atom_id_2 + 1, atom_id_3 + 1))
    out_file.write('\n')
    out_file.close()

i = 0
for y in range(1):
    for z in range(1):
        n = np.array([1, 0, 0])#np.random.randn(3)
        add_nanoparticle(np.array([0,0,0]), 20, 20, 24, i)
        #build_chain(i, 24, np.array([0,0,0]), n / np.linalg.norm(n, axis=0))
        i += 1

#add_nanoparticles(10, i)
max_x *= 10
max_y *= 10
max_z *= 10
output_polymer()
