# LAMMPS Experiments
```setup_env.sh``` scripts sets up lammps install via conda, creating a new 'lammps' environment, which can be used using:
```conda activate lammps```


VMD can be used to visualize results: https://www.ks.uiuc.edu/Research/vmd/alpha/

tut01 is from: https://lammpstutorials.github.io/tutorials/tutorial01.html
Output file is dump.lammpstrj, visualize using:
```
vmd dump.lammpstrj
```

## Building for Cori

```bash
module load cmake
module load intel
module load openmpi
git clone -b release https://github.com/lammps/lammps.git mylammps
cd mylammps
mkdir build
cd build
wget -O voro_build-prefix/src/voro++-0.4.6.tar.gz https://download.lammps.org/thirdparty/voro++-0.4.6.tar.gz
wget -O Eigen3_build-prefix/src/eigen-3.4.0.tar.gz https://download.lammps.org/thirdparty/eigen-3.4.0.tar.gz
cmake -C ../cmake/presets/most.cmake ../cmake -D BUILD_MPI=yes -D LAMMPS_MACHINE=custom -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_Fortran_COMPILER=mpifort -DPKG_BROWNIAN=yes
make -j32
make install
```

Path to executable:
```
~/.local/bin/lmp_custom
```