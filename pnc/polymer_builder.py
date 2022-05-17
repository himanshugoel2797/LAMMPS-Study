# Copyright (c) 2022 Himanshu Goel
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

import os
import sys

chain_length = 86 # corresponds to 9kg/mol
one_d = 5.29 # Angstrom
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
def build_monomer(chain_idx, monomer_idx, chain_pos_y=0, chain_pos_z=0):
    global max_x
    global max_y
    global max_z
    x_pos = monomer_idx
    b_id = len(atoms)
    h_id = len(atoms) + 1
    atoms.append((b_id, chain_idx, 1, x_pos, chain_pos_y, chain_pos_z))
    atoms.append((h_id, chain_idx, 2, x_pos, chain_pos_y + H_dist, chain_pos_z)) #h bead is 90 degrees from b bead

    if max_x < x_pos:
        max_x = x_pos
    if max_y < chain_pos_y:
        max_y = chain_pos_y
    if max_z < chain_pos_z:
        max_z = chain_pos_z

    #Bonds
    if monomer_idx > 0: #bond to previous monomer
        bonds.append((len(bonds), 1, b_id - 2, b_id))
    bonds.append((len(bonds), 2, b_id, h_id)) #bond to h bead

    #Angles
    if monomer_idx > 0: #angle to previous monomer
        angles.append((len(angles), 1, b_id - 2, b_id, h_id))

def add_nanoparticles(cnt, chain_idx):
    for i in range(cnt):
        x_pos = i * (box_side / cnt)
        y_pos = 0
        z_pos = 0
        atoms.append((len(atoms), chain_idx, 3, x_pos, y_pos, z_pos))

def build_chain(chain_idx, chain_pos_y=0, chain_pos_z=0):
    for i in range(chain_length):
        build_monomer(chain_idx, i, chain_pos_y, chain_pos_z)

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
    out_file.write('{} bond types\n'.format(2))
    out_file.write('{} angle types\n'.format(1))
    out_file.write('\n')
    out_file.write('{} {} xlo xhi\n'.format(0, (max_x + 1) * one_d))
    out_file.write('{} {} ylo yhi\n'.format(0, (max_y + 1) * one_d))
    out_file.write('{} {} zlo zhi\n'.format(0, (max_z + 1) * one_d))
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
    for atom_id, chain_idx, atom_type, x_pos, y_pos, z_pos in atoms:
        out_file.write('{} {} {} {} {} {}\n'.format(atom_id + 1, chain_idx + 1, atom_type, (x_pos + 0.5) * one_d, (y_pos + 0.5) * one_d, (z_pos + 0.5) * one_d))
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
for y in range(15):
    for z in range(15):
        build_chain(i, chain_pos_y = y * 2, chain_pos_z = z * 2)
        i += 1

#add_nanoparticles(10, i)
max_x *= 1.2
max_y *= 1.2
max_z *= 1.2
output_polymer()
