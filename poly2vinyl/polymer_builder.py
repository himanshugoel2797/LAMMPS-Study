# Copyright (c) 2022 Himanshu Goel
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

import os
import sys

chain_length = 36
one_d = 5.29 # Angstrom
H_d = 0.3 #d units
H_dist = 0.37 #d units
H_mass = 17.007 #g/mol for OH group
B_mass = 105.14 - H_mass #g/mol for the r est of the monomer

atoms = []
bonds = []
angles = []
dihedrals = []

def build_monomer(chain_idx):
    x_pos = len(atoms)
    b_id = len(atoms)
    h_id = len(atoms) + 1
    atoms.append((b_id, chain_idx, 1, x_pos, 0, 0))
    atoms.append((h_id, chain_idx, 2, x_pos, H_dist, 0)) #h bead is 90 degrees from b bead

    #Bonds
    if b_id - 2 >= 0: #bond to previous monomer
        bonds.append((len(bonds), 1, b_id - 2, b_id))
    bonds.append((len(bonds), 2, b_id, h_id)) #bond to h bead

    #Angles
    if b_id - 2 >= 0: #angle to previous monomer
        angles.append((len(angles), 1, b_id - 2, b_id, h_id))

def build_chain(chain_idx):
    for i in range(chain_length):
        build_monomer(chain_idx)

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
    out_file.write('{} atom types\n'.format(2))
    out_file.write('{} bond types\n'.format(2))
    out_file.write('{} angle types\n'.format(1))
    out_file.write('\n')
    out_file.write('{} {} xlo xhi\n'.format(0, (chain_length + 1) * one_d))
    out_file.write('{} {} ylo yhi\n'.format(0, 2 * (H_dist + H_d) * one_d))
    out_file.write('{} {} zlo zhi\n'.format(0, 1 * one_d))
    out_file.write('\n')
    out_file.write('\n')
    out_file.write('Masses\n')
    out_file.write('\n')
    out_file.write('1 {}\n'.format(B_mass))
    out_file.write('2 {}\n'.format(H_mass))
    out_file.write('\n')
    out_file.write('\n')
    out_file.write('Atoms\n')
    out_file.write('\n')
    for atom_id, chain_idx, atom_type, x_pos, y_pos, z_pos in atoms:
        out_file.write('{} {} {} {} {} {}\n'.format(atom_id, chain_idx, atom_type, x_pos * one_d, y_pos * one_d, z_pos * one_d))
    out_file.write('\n')
    out_file.write('\n')
    out_file.write('Bonds\n')
    out_file.write('\n')
    for bond_id, bond_type, atom_id_1, atom_id_2 in bonds:
        out_file.write('{} {} {} {}\n'.format(bond_id, bond_type, atom_id_1, atom_id_2))
    out_file.write('\n')
    out_file.write('\n')
    out_file.write('Angles\n')
    out_file.write('\n')
    for angle_id, angle_type, atom_id_1, atom_id_2, atom_id_3 in angles:
        out_file.write('{} {} {} {} {}\n'.format(angle_id, angle_type, atom_id_1, atom_id_2, atom_id_3))
    out_file.write('\n')
    out_file.close()

build_chain(0)
output_polymer()