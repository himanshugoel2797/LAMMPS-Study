# Generates version of script with given concentration and associated slurm script

import os
import sys
import numpy as np

src_file = 'bd_v3.lammps'
dest_dir = 'concent_{}'
dest_file_path = '{}/bd.lammps'

vol = 400000. * 400000. * 100000.
particle_diameter = 2500

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

concents = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
for idx, i in enumerate(concents):
    particle_cnt = generate_updated_file(i)
    print('{}/{} [{} particles]'.format(idx, len(concents), particle_cnt))
    if idx != 0 and idx % 2 == 0:
        # Wait for 30 minutes
        print('Waiting for 30 minutes')
        os.system('sleep 30m')