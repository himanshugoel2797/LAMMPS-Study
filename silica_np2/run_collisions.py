# Generates version of script with given concentration and associated slurm script

import os
import sys
import numpy as np

src_file = 'bd_v3_collisions.lammps'
dest_dir = 'concent_col_{}'
dest_file_path = '{}/bd.lammps'

vol = 400000. * 400000. * 400000.
particle_diameter = 2500

def parse_bound(line):
    line = line.split(' ', 1)
    return float(line[0]), float(line[1])

def parse_atom_def(line):
    line = line.split(' ')
    return int(line[0]), int(line[1]), float(line[2]), float(line[3]), float(line[4])

def parse_timestep(f):
    if len(f.readline()) == 0: #Check for new timestamp
        return None
    timestep = int(f.readline())

    f.readline()  # Skip num of atoms
    atomCount = int(f.readline())

    f.readline()  # Skip box bounds
    xlo, xhi = parse_bound(f.readline())
    ylo, yhi = parse_bound(f.readline())
    zlo, zhi = parse_bound(f.readline())

    f.readline() # Skip atom definitions
    atomDefs = []
    for i in range(atomCount):
        atomDefs.append(parse_atom_def(f.readline()))

    return timestep, atomCount, xlo, xhi, ylo, yhi, zlo, zhi, atomDefs

def parse_file(filename):
    f = open(filename, 'r')
    timesteps = []
    while True:
        timestep = parse_timestep(f)
        if timestep is None:
            break
        timesteps.append(timestep)
    f.close()
    return timesteps

def convert_file(filename, output_pattern, atom_radii, sc = 5e-7):
    timesteps = parse_file(filename)
    for timestep in timesteps:
        timestep_filename = output_pattern.format(timestep[0])
        f = open(timestep_filename, 'w+')

        diff_x = timestep[3] - timestep[2]
        add_x = timestep[2]

        diff_y = timestep[5] - timestep[4]
        add_y = timestep[4]

        diff_z = timestep[7] - timestep[6]
        add_z = timestep[6]

        for atom in timestep[8]:
            f.write('{}\t{}\t{}\tS\t{}\n'.format((atom[2] * diff_x + add_x) * sc, (atom[3] * diff_y + add_y) * sc, (atom[4] * diff_z + add_z) * sc, atom_radii[atom[1] - 1]))
        f.close()

def generate_updated_file(concent=0.02):
    # Calculate particle count for given concentration
    particle_count = int(vol * (concent / (4./3. * np.pi * (particle_diameter*0.5)**3)))
    # Read src_file and replace 'REPLACE_ME' with particle count
    with open(src_file, 'r') as f:
        content = f.read()
    content = content.replace('REPLACE_ME', str(particle_count)).replace('CONCENT', str(concent))
    # Write to dest_file
    dest_dir_l = dest_dir.format(concent)
    if not os.path.exists(dest_dir_l):
        os.makedirs(dest_dir_l)
    dest_file = dest_file_path.format(dest_dir_l)
    with open(dest_file, 'w') as f:
        f.write(content)

    #Copy slurm script
    slurm_file = '{}/cori.slurm'.format(dest_dir_l)
    with open('cori.slurm', 'r') as f:
        content = f.read()
    content = content.replace('REPLACE_ME', str(particle_count))
    with open(slurm_file, 'w') as f:
        f.write(content)
    # change working directory to dest_dir_l
    os.chdir(dest_dir_l)
    os.system('sbatch cori.slurm')
    os.chdir('..')
    return particle_count

concents = [0.02, 0.2]
for idx, i in enumerate(concents):
    particle_cnt = generate_updated_file(i)

    print('{}/{} [{} particles]'.format(idx, len(concents), particle_cnt))
#    if idx != 0 and idx % 2 == 0:
        # Wait for 30 minutes
#        print('Waiting for 30 minutes')
#        os.system('sleep 30m')

#for idx, concent in enumerate(concents):
#    out_dir = 'converted_{}'.format(concent)
#    #Make directory for converted files
#    if not os.path.exists(out_dir):
#        os.makedirs(out_dir)
#
#    convert_file('concent_{}/dump_{}.lammpstrj'.format(concent, concent), out_dir + '/smp_{}.dat', [125.e-9], 20.e-6)
