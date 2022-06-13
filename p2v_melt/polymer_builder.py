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
chain_surface_density = 0.2
box_side = 500 #529nm^3
NP_count = 8
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
nanoparticles = []

min_pos = np.zeros(3)
max_pos = np.zeros(3)

def sample_cube(samples):
    side = int(math.ceil(samples ** (1./3.)))

    points = []
    for z in range(side):
        for y in range(side):
            for x in range(side):
                if len(points) < samples:
                    points.append(np.array([x/side, y/side, z/side]))
                else:
                    return points
    return points

def sample_spherical(samples):
    if samples == 1:
        return np.array([[0, 0, 0]])

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
    global min_pos
    global max_pos

    b_id = len(atoms)
    h_id = len(atoms) + 1

    #generate random 3d vector using numpy
    b_pos = origin + direction * 2 * B_rad * monomer_idx
    h_pos = b_pos + perpendicular_vector(direction) * H_dist

    for i in range(3):
        min_pos[i] = min([min_pos[i], b_pos[i], h_pos[i]])
        max_pos[i] = max([max_pos[i], b_pos[i], h_pos[i]])

    if monomer_idx == 0:
        atoms.append((b_id, chain_idx, 4, b_pos))
    else:
        atoms.append((b_id, 0, 1, b_pos))
    atoms.append((h_id, 0, 2, h_pos)) #h bead is 90 degrees from b bead
        
    #Bonds
    if monomer_idx > 0: #bond to previous monomer
        bonds.append((len(bonds), 1, b_id - 2, b_id))
    bonds.append((len(bonds), 2, b_id, h_id)) #bond to h bead

    #Angles
    if monomer_idx > 0: #angle to previous monomer
        angles.append((len(angles), 1, h_id - 2, b_id - 2, b_id))

def build_freemonomer(monomer_idx, origin, direction):
    global min_pos
    global max_pos

    b_id = len(atoms)
    h_id = len(atoms) + 1

    #generate random 3d vector using numpy
    b_pos = origin + direction * 2 * B_rad
    h_pos = b_pos + perpendicular_vector(direction) * H_dist

    for i in range(3):
        min_pos[i] = min([min_pos[i], b_pos[i], h_pos[i]])
        max_pos[i] = max([max_pos[i], b_pos[i], h_pos[i]])

    atoms.append((b_id, 0, 1, b_pos))
    atoms.append((h_id, 0, 2, h_pos)) #h bead is 90 degrees from b bead
        
    #Bonds
    if monomer_idx > 0: #bond to previous monomer
        bonds.append((len(bonds), 1, b_id - 2, b_id))
    bonds.append((len(bonds), 2, b_id, h_id)) #bond to h bead

    #Angles
    if monomer_idx > 0: #angle to previous monomer
        angles.append((len(angles), 1, h_id - 2, b_id - 2, b_id))
    return b_pos


def build_boundchain(chain_idx, chain_len, origin, direction):
    for i in range(chain_len):
        build_monomer(chain_idx, i, origin, direction)

def build_randomchain(chain_idx, chain_len, origin):
    for i in range(chain_len):
        direction = np.random.rand(3)
        direction /= np.linalg.norm(direction)
        origin = build_freemonomer(chain_idx, i, origin, direction)

def add_nanoparticle(pos, rad, chain_cnt, chain_len, chain_idx):
    np_idx = len(atoms)
    atoms.append((np_idx, chain_idx, 3, pos))
    
    atom_idxs = [np_idx]

    origins = sample_spherical(chain_cnt)
    for i in range(chain_cnt):
        atom_idxs.append(len(atoms))
        #bonds.append((len(bonds), 3, np_idx, len(atoms)))
        build_boundchain(chain_idx, np.random.randint(12, chain_len+1), pos + origins[i] * (rad + B_rad), origins[i])

    nanoparticles.append((len(nanoparticles), chain_idx, atom_idxs))


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
    out_file.write('{} {} xlo xhi\n'.format(min_pos[0] - 5, max_pos[0] + 5))
    out_file.write('{} {} ylo yhi\n'.format(min_pos[1] - 5, max_pos[1] + 5))
    out_file.write('{} {} zlo zhi\n'.format(min_pos[2] - 5, max_pos[2] + 5))
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
        out_file.write('{} {} {} {} {} {}\n'.format(atom_id + 1, chain_idx, atom_type, pos[0] * one_d, pos[1] * one_d, pos[2] * one_d))
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

def output_def(filename='polymer_def_builder.lammps'):
    out_file = open(filename, 'w')
    out_file.write('# Parameters and groups for poly(2-vinyl) chain\n')
    out_file.write('\n')
    out_file.write('variable np_rad equal {}\n'.format(NP_rad))

    molecule_ids = ""
    for _, id, _ in nanoparticles:
        molecule_ids += str(id) + " "
    out_file.write('group nanoparticles molecule {}\n'.format(molecule_ids))
    out_file.write('group polymers molecule 0\n')
    out_file.close()
    
pos_set = sample_cube(NP_count + 500)
min_dist = math.inf
if NP_count > 1:
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
    add_nanoparticle(pos, NP_rad, chain_surface_count, chain_length, i + 1)

output_polymer()
output_def()
