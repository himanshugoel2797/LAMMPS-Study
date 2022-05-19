NVARCH=`uname -s`_`uname -m`; export NVARCH
NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
NVVERSION='22.3'; export NVVERSION
MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/$NVVERSION/compilers/man; export MANPATH
PATH=$NVCOMPILERS/$NVARCH/$NVVERSION/compilers/bin:$PATH; export PATH
LD_LIBRARY_PATH=$NVCOMPILERS/$NVARCH/$NVVERSION/math_libs/lib64:$NVCOMPILERS/$NVARCH/$NVVERSION/cuda/lib64:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH