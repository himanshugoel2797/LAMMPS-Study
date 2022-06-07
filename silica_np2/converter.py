import os
import readline

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

concent = 0.02
out_dir = 'converted_{}'.format(concent)
#Make directory for converted files
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

convert_file('concent_{}/dump_{}.lammpstrj'.format(concent, concent), out_dir + '/smp_{}.dat', [125.e-9], 20.e-6)
