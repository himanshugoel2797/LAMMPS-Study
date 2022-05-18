# Copyright (c) 2022 Himanshu Goel
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

#Generates a polymer nanocomposite of poly(2-vinylpyridine) with a nanoparticle of SiO2

import os
import sys
import math
import numpy as np

chain_length = 86 # corresponds to 9kg/mol
chain_surface_density = 0.05
box_side = 500 #529nm^3
NP_count = 10
NP_rad_nm = 9.1 #nm

one_d = 1#5.29 # Angstrom
H_rad = 0.15 #d units
B_rad = 0.5 #d units
H_dist = 0.37 #d units
H_mass = 0.001
B_mass = 1 - H_mass #B-beads are treated as center of mass of the monomer

NP_rad = NP_rad_nm * 10 / 5.29 #d units
NP_volume = (4/3) * np.pi * NP_rad**3
NP_surface_area = np.pi * NP_rad**2
NP_density = 0.01949
NP_mass = NP_density * NP_volume
chain_surface_count = int(NP_surface_area * chain_surface_density)

print ("Nanoparticle count:", NP_count)
print ("Chains per nanoparticle:", chain_surface_count)
print ("Nanoparticle size in d:", NP_rad)
res = input("Generate? (y/n)")
if res.capitalize() != "Y":
    sys.exit()

atoms = []
bonds = []
angles = []
dihedrals = []

def sample_spherical(samples=1000):

    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append(np.array([x, y, z]))

    return np.array(points)

def perpendicular_vector(v):
    if v[1] == 0 and v[2] == 0:
        if v[0] == 0:
            raise ValueError('zero vector')
        else:
            return np.cross(v, [0, 1, 0])
    return np.cross(v, [1, 0, 0])

def build_monomer(chain_idx, monomer_idx, origin, direction):
    b_id = len(atoms)
    h_id = len(atoms) + 1

    #generate random 3d vector using numpy
    b_pos = origin + direction * 2 * B_rad * monomer_idx
    h_pos = b_pos + perpendicular_vector(direction) * H_dist

    if monomer_idx == 0:
        atoms.append((b_id, chain_idx, 4, b_pos))
    else:
        atoms.append((b_id, chain_idx, 1, b_pos))
    atoms.append((h_id, chain_idx, 2, h_pos)) #h bead is 90 degrees from b bead
        
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
    
    origins = sample_spherical(chain_cnt)
    for i in range(chain_cnt):
        bonds.append((len(bonds), 3, np_idx, len(atoms)))
        build_chain(chain_idx, np.random.randint(2, chain_len+1), pos + origins[i] * (rad + B_rad), origins[i])


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
    out_file.write('{} atom types\n'.format(4))
    out_file.write('{} bond types\n'.format(3))
    out_file.write('{} angle types\n'.format(1))
    out_file.write('\n')
    out_file.write('{} {} xlo xhi\n'.format(-box_side, box_side))
    out_file.write('{} {} ylo yhi\n'.format(-box_side, box_side))
    out_file.write('{} {} zlo zhi\n'.format(-box_side, box_side))
    out_file.write('\n')
    out_file.write('\n')
    out_file.write('Masses\n')
    out_file.write('\n')
    out_file.write('1 {}\n'.format(B_mass))
    out_file.write('2 {}\n'.format(H_mass))
    out_file.write('3 {}\n'.format(NP_mass))
    out_file.write('4 {}\n'.format(B_mass))
    out_file.write('\n')
    out_file.write('\n')
    out_file.write('Atoms\n')
    out_file.write('\n')
    for atom_id, chain_idx, atom_type, pos in atoms:
        out_file.write('{} {} {} {} {} {}\n'.format(atom_id + 1, chain_idx + 1, atom_type, pos[0] * one_d, pos[1] * one_d, pos[2] * one_d))
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

pos_set = sample_spherical(NP_count)
min_dist = math.inf
for i in range(NP_count):
    for j in range(i+1, NP_count):
        dist = np.linalg.norm(pos_set[i] - pos_set[j])
        if dist < min_dist:
            min_dist = dist
#scale pos_set such that min_dist is NP_rad
for i in range(NP_count):
    pos_set[i] = pos_set[i] * 2 * (NP_rad + B_rad * 2 * chain_length) / min_dist
box_side = np.max(pos_set) + (NP_rad + B_rad * 2 * chain_length)

for i in range(NP_count):
    pos = pos_set[i]
    add_nanoparticle(pos, NP_rad, chain_surface_count, chain_length, i)

output_polymer()
