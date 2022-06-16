/**
 * Copyright (c) 2022 Himanshu Goel
 * 
 * This software is released under the MIT License.
 * https://opensource.org/licenses/MIT
 */

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <tuple>

const int grid_side = 100;
const int grid_size = grid_side * grid_side * grid_side;
const int sparsity = 1;
const float surface_chain_prob = 0.4f;
const bool surface_chains = true;

const float grid_scale = 1;//1.e-9; // each voxel is 0.5nm

const int NP_COUNT = 50;
const float NP_RADIUS_UNSCALED = 4;//9.0e-9;
const float NP_RADIUS = NP_RADIUS_UNSCALED / grid_scale;
const char NP_TYPE = 1;
const double NP_MASS = 64;//2650 * (4.0 / 3.0) * M_PI * pow(NP_RADIUS_UNSCALED, 3);

const int BEAD_LENGTH_MAX = 24;
const float BEAD_RADIUS = 0.5f;//0.5e-9 / grid_scale;
const char BEAD_TYPE_A = 2;
const char BEAD_TYPE_B = 3;
const double BEAD_MASS = 1;//(0.10514 * 6.02214e-23);

int grid[grid_side][grid_side][grid_side];

#define grid_get(x, y, z) grid[x][y][z]
#define grid_set_zr(x, y, z, v) grid[x][y][z] = v

struct particle_def_t {
    public:
        int x, y, z;
        int chain_idx;
        char type;
    
        particle_def_t(int x, int y, int z, char type, int chain_idx = 0) {
            this->x = x;
            this->y = y;
            this->z = z;
            this->type = type;
            this->chain_idx = chain_idx;
        }
        particle_def_t(){}
};
struct bond_def_t {
    public:
        int a, b;
        char type;
    
        bond_def_t(int a, int b, char type) {
            this->a = a;
            this->b = b;
            this->type = type;
        }
        bond_def_t(){}
};
struct angle_def_t {
    public:
        int a, b, c;
    
        angle_def_t(int a, int b, int c) {
            this->a = a;
            this->b = b;
            this->c = c;
        }
        angle_def_t(){}
};

int add_np(int x, int y, int z, std::vector<particle_def_t> &particles, std::vector<bond_def_t> &bonds, std::vector<angle_def_t> &angles) {
    if (x < 0 || x >= grid_side || y < 0 || y >= grid_side || z < 0 || z >= grid_side) {
        return 0;
    }

    int s = 2;
    int t = 0;

    // Make sure volume is empty
    for (int i = -NP_RADIUS - s; i <= NP_RADIUS + s; i++) {
        for (int j = -NP_RADIUS - s; j <= NP_RADIUS + s; j++) {
            for (int k = -NP_RADIUS - s; k <= NP_RADIUS + s; k++) {
                if (i * i + j * j + k * k <= ceil((NP_RADIUS + t) * (NP_RADIUS + t))) {
                    if (x + i < 0 || x + i >= grid_side || y + j < 0 || y + j >= grid_side || z + k < 0 || z + k >= grid_side) {
                        printf ("Boundary rejection %d %d %d.\n", x+i, y+j, z+k);
                        return 0;
                    }
                    if (grid_get(x + i, y + j, z + k)) {
                        printf ("Overlap rejection %d %d %d %d %d %d %d.\n", x, y, z, i, j, k, grid_get(x + i, y + j, z + k));
                        return 0;
                    }
                }
            }
        }
    }

    int p_idx = particles.size();

    // Fill volume
    for (int i = -NP_RADIUS - s; i <= NP_RADIUS + s; i++) {
        for (int j = -NP_RADIUS - s; j <= NP_RADIUS + s; j++) {
            for (int k = -NP_RADIUS - s; k <= NP_RADIUS + s; k++) {
                if (i * i + j * j + k * k <= ceil((NP_RADIUS + t) * (NP_RADIUS + t))) {
                    grid_set_zr(x + i, y + j, z + k, p_idx);
                }
            }
        }
    }

    // Add particle
    particle_def_t np(x, y, z, NP_TYPE);
    particles.push_back(np);
    
    return 1;
}

int bond_surface_chain(int np_idx, std::vector<particle_def_t> &particles, std::vector<bond_def_t> &bonds, std::vector<angle_def_t> &angles, std::vector<int> &bonded_chains) {
    int x = particles[np_idx].x;
    int y = particles[np_idx].y;
    int z = particles[np_idx].z;

    int s = 2;
    int t = 1;

    std::vector<int> bonded_chain_idxs;

    // Find surface chains to bond to
    for (int i = -NP_RADIUS - s; i <= NP_RADIUS + s; i++) {
        for (int j = -NP_RADIUS - s; j <= NP_RADIUS + s; j++) {
            for (int k = -NP_RADIUS - s; k <= NP_RADIUS + s; k++) {
                if (i * i + j * j + k * k <= ceil((NP_RADIUS + t) * (NP_RADIUS + t))) {
                    if (x + i < 0 || x + i >= grid_side || y + j < 0 || y + j >= grid_side || z + k < 0 || z + k >= grid_side)
                        continue;

                    if (grid_get(x + i, y + j, z + k)) {
                        int tgt_idx = grid_get(x + i, y + j, z + k);
                        if (particles[tgt_idx].type == BEAD_TYPE_A && rand() < RAND_MAX * surface_chain_prob) {
                            bool already_bonded = std::find(bonded_chain_idxs.begin(), bonded_chain_idxs.end(), particles[tgt_idx].chain_idx) != bonded_chain_idxs.end();
                            if (!already_bonded)
                                already_bonded = std::find(bonded_chains.begin(), bonded_chains.end(), particles[tgt_idx].chain_idx) != bonded_chains.end();
                            if (already_bonded)
                                continue;
                            bond_def_t bond(np_idx, tgt_idx, 1);
                            bonds.push_back(bond);
                            bonded_chain_idxs.push_back(particles[tgt_idx].chain_idx);
                            bonded_chains.push_back(particles[tgt_idx].chain_idx);
                        }
                    }
                }
            }
        }
    }

    return 1;
}

int add_free_polymer(int x, int y, int z, int chain_idx, std::vector<particle_def_t> &particles, std::vector<bond_def_t> &bonds, std::vector<angle_def_t> &angles) {
    if (x < 0 || x >= grid_side || y < 0 || y >= grid_side || z < 0 || z >= grid_side) {
        return 0;
    }

    // Place bead
    particles.push_back(particle_def_t(x, y, z, BEAD_TYPE_A, chain_idx));
    int prev_bead_index = particles.size() - 1;
    grid_set_zr(x, y, z, prev_bead_index);

    // Place bonds
    for (int i = 0; i < BEAD_LENGTH_MAX; i++) {
        std::vector<std::tuple<int, int, int>> available_dirs;
        for (int j = 0; j < 3; j++) {
            int x_dir = 0, y_dir = 0, z_dir = 0;
            for (int j_sign = -1; j_sign <= 1; j_sign+=2) {
                x_dir = j == 0 ? j_sign : 0;
                y_dir = j == 1 ? j_sign : 0;
                z_dir = j == 2 ? j_sign : 0;
                
                if (x + x_dir < 0 || x + x_dir >= grid_side || y + y_dir < 0 || y + y_dir >= grid_side || z + z_dir < 0 || z + z_dir >= grid_side) {
                    continue;
                }

                if (!grid_get(x + x_dir, y + y_dir, z + z_dir)) {
                    available_dirs.push_back(std::make_tuple(x_dir, y_dir, z_dir));
                }
            }
        }

        if (available_dirs.size() == 0) {
            //printf ("Dead end termination %d at length %d.\n", chain_idx, i);
            break;
        }

        int x_dir, y_dir, z_dir;
        std::tie(x_dir, y_dir, z_dir) = available_dirs[rand() % available_dirs.size()];
        int new_x = x + x_dir;
        int new_y = y + y_dir;
        int new_z = z + z_dir;
        
        particles.push_back(particle_def_t(new_x, new_y, new_z, BEAD_TYPE_A, chain_idx));
        bonds.push_back(bond_def_t(prev_bead_index, particles.size() - 1, 0));
        prev_bead_index = particles.size() - 1;
        grid_set_zr(new_x, new_y, new_z, prev_bead_index);

        x = new_x;
        y = new_y;
        z = new_z;
    }
    return 1;
}

int main(){
    for (int i = 0; i < grid_side; i++)
        for (int j = 0; j < grid_side; j++)
            for (int k = 0; k < grid_side; k++)
                grid_set_zr(i, j, k, 0);
    std::vector<particle_def_t> particles;
    std::vector<bond_def_t> bonds;
    std::vector<bond_def_t> surface_bonds;
    std::vector<angle_def_t> angles;
    std::vector<int> bonded_chains;

    // Add dummy particle
    particles.push_back(particle_def_t(0, 0, 0, 0));

    // Set fixed seed
    srand(0);

    //Place the nanoparticles in the grid with their surface polymers
    for (int i = 0; i < NP_COUNT; i++) {
        while (!add_np(rand() % grid_side, rand() % grid_side, rand() % grid_side, particles, bonds, angles))
            printf("Failed to add nanoparticle %d, Trying again.\n", i);
    }

    //Fill the remaining space with free polymers
    //Find the first empty voxel and start a polymer there
    int chain_idx = 1;
    for (int i = 0; i < grid_side; i+=sparsity) {
        for (int j = 0; j < grid_side; j+=sparsity) {
            for (int k = 0; k < grid_side; k+=sparsity) {
                if (!grid_get(i, j, k)) {
                    //if (chain_idx < 300){
                    if (add_free_polymer(i, j, k, chain_idx, particles, bonds, angles))
                        chain_idx++;
                    if (chain_idx % 10000 == 0)
                        printf("Generated chain #%d\n", chain_idx);
                    //}
                }
            }
        }
    }

    if (surface_chains) {
        for (int i = 0; i < NP_COUNT; i++) {
            bond_surface_chain(i+1, particles, surface_bonds, angles, bonded_chains);
        }
        bonds.insert(bonds.begin(), surface_bonds.begin(), surface_bonds.end());
    }

    //Write the data to a LAMMPS input file
    FILE *fp = fopen("pnc_def.lammps", "w");
    fprintf(fp, "# LAMMPS data file for PNC\n");
    fprintf(fp, "%d atoms\n", particles.size() - 1);
    if (bonds.size() > 0) fprintf(fp, "%d bonds\n", bonds.size());
    //fprintf(fp, "%d angles\n", angles.size());
    fprintf(fp, "\n\n");
    fprintf(fp, "%d atom types\n", surface_chains ? 3 : 2);
    if (bonds.size() > 0) fprintf(fp, "%d bond types\n", 2);
    //fprintf(fp, "%d angle types\n", ANGLE_TYPE + 1);
    fprintf(fp, "\n\n");
    fprintf(fp, "%g %g xlo xhi\n", -0.01 * grid_side * grid_scale, (grid_side * 1.01) * grid_scale);
    fprintf(fp, "%g %g ylo yhi\n", -0.01 * grid_side * grid_scale, (grid_side * 1.01) * grid_scale);
    fprintf(fp, "%g %g zlo zhi\n", -0.01 * grid_side * grid_scale, (grid_side * 1.01) * grid_scale);
    fprintf(fp, "\n\n");
    fprintf(fp, "Masses\n");
    fprintf(fp, "\n");
    fprintf(fp, "%d %.20g\n", NP_TYPE, NP_MASS);
    fprintf(fp, "%d %.20g\n", BEAD_TYPE_A, BEAD_MASS);
    if (surface_chains)
        fprintf(fp, "%d %.20g\n", BEAD_TYPE_B, BEAD_MASS);
    fprintf(fp, "\n\n");
    fprintf(fp, "Atoms\n");
    fprintf(fp, "\n");
    for (int i = 1; i < particles.size(); i++) {
        bool chain_bonded = std::find(bonded_chains.begin(), bonded_chains.end(), particles[i].chain_idx) != bonded_chains.end();
        if (chain_bonded)
            fprintf(fp, "%d %d %d %g %g %g\n", i, particles[i].chain_idx, BEAD_TYPE_B, particles[i].x * grid_scale, particles[i].y * grid_scale, particles[i].z * grid_scale);
        else
            fprintf(fp, "%d %d %d %g %g %g\n", i, particles[i].chain_idx, particles[i].type, particles[i].x * grid_scale, particles[i].y * grid_scale, particles[i].z * grid_scale);
    }
    fprintf(fp, "\n\n");
    if (bonds.size() > 0){
        fprintf(fp, "Bonds\n");
        fprintf(fp, "\n");
        for (int i = 0; i < bonds.size(); i++) {
            fprintf(fp, "%d %d %d %d\n", i + 1, bonds[i].type + 1, bonds[i].a, bonds[i].b);
        }
        fprintf(fp, "\n\n");
    }
    //fprintf(fp, "Angles\n");
    //fprintf(fp, "\n");
    //for (int i = 0; i < angles.size(); i++) {
    //    fprintf(fp, "%d %d %d %d %d\n", i + 1, 1, angles[i].a + 1, angles[i].b + 1, angles[i].c + 1);
    //}
    //fprintf(fp, "\n");
    fclose(fp);
}