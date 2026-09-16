# ChihayaAnon: GPU accelerated ShengBTE

Install(the same as OpenACC version VASP):
1. Download and install intel oneapi-toolkit
2. Download and install NVIDIA HPC-SDK
3. Download and install SPGLIB
4. Change the spglib path in arch.make to your spglib path
5. Activate oneapi environment first(To use its MKL)
6. Activate NV HPC-SDK environment (check by entering which mpifort or which mpirun)
7. Change the -gpu=cc70 in arch.make to your GPU compute compability(cc70 is V100 cc89 is RTX4090)
8. make

## Usage

The force-constant input format is the same as that used by ShengBTE.

This version introduces several additional options in the configuration file.

### Flags

- `itergpu`: Default is `True`. Performs the iterative calculation on the GPU. If many CPU cores are available, setting this option to `False` may provide better performance.
- `DIIS`: Default is `True`. Enables damped Pulay mixing.

### Parameters

- `mmax`: Default is `6`. Specifies the maximum number of history arrays stored by DIIS. It generally does not need to be changed; if the calculation fails to converge, try increasing this value.
- `beta_pulay`: Default is `0.2`. Specifies the damping coefficient for DIIS. It generally does not need to be changed.
- `alpha`: Specifies the linear-mixing coefficient used when DIIS is disabled. The default value is `0.5`, which mixes the old and newly calculated results equally.

### Example input

```fortran
&allocations
    nelements=4
    natoms=24
    ngrid(:)=12 12 12
&end

&crystal
    lfactor=0.1000000000
    lattvec(:,1)=10.393430565951 0.000000037792 0.000000577490
    lattvec(:,2)=-0.000000019342 3.916346936511 0.000000387192
    lattvec(:,3)=0.000000842181 0.000001460396 14.488396423858
    elements="Lu" "S" "Tl" "Zn"
    types=1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 4 4 4 4
    positions(:,1)=0.486382859227 0.749999981887 0.247967790494
    positions(:,2)=0.986382860167 0.750000015759 0.252032195316
    positions(:,3)=0.013617141588 0.249999984442 0.747967802100
    positions(:,4)=0.513617138058 0.250000018932 0.752032208294
    positions(:,5)=0.425511787296 0.250000006893 0.128912170989
    positions(:,6)=0.045874636230 0.250000034428 0.131137178890
    positions(:,7)=0.737000340219 0.750000204543 0.189011055950
    positions(:,8)=0.237000337800 0.749999792384 0.310988933728
    positions(:,9)=0.545874617853 0.249999964018 0.368862806070
    positions(:,10)=0.925511767261 0.249999991099 0.371087817244
    positions(:,11)=0.074488234460 0.750000007307 0.628912183461
    positions(:,12)=0.454125381866 0.750000036587 0.631137196513
    positions(:,13)=0.762999661978 0.250000207768 0.689011066793
    positions(:,14)=0.262999659356 0.249999794859 0.810988943964
    positions(:,15)=0.954125362847 0.749999966372 0.868862819683
    positions(:,16)=0.574488214661 0.749999995521 0.871087828692
    positions(:,17)=0.765120213580 0.250000209695 0.001638641996
    positions(:,18)=0.265120257666 0.249999790096 0.498361362453
    positions(:,19)=0.734879742092 0.750000208475 0.501638637626
    positions(:,20)=0.234879786355 0.749999789697 0.998361358034
    positions(:,21)=0.237026852391 0.249999687449 0.213647470077
    positions(:,22)=0.737026870243 0.250000309309 0.286352518106
    positions(:,23)=0.262973131233 0.749999689410 0.713647480796
    positions(:,24)=0.762973145574 0.750000313070 0.786352532733
    scell(:)=2 2 2
&end

&parameters
    T=300
    scalebroad=0.1
&end

&flags
    nonanalytic=.FALSE.
    nanowires=.FALSE.
    convergence=.TRUE.
    itergpu=.TRUE.
    DIIS=.TRUE.
&end

```

### Running 

```bash
#activate NV-HPC environment first 
export NVCOMPILERS=/home14/liangqm/QUANTUM_ESPRESSO/hpc_sdk/Linux_x86_64
export NVHPC_VERSION=24.5

export PATH=$NVCOMPILERS/$NVHPC_VERSION/compilers/bin:$PATH
export PATH=$NVCOMPILERS/$NVHPC_VERSION/comm_libs/mpi/bin:$PATH

export MANPATH=$MANPATH:$NVCOMPILERS/$NVHPC_VERSION/compilers/man
export MANPATH=$MANPATH:$NVCOMPILERS/$NVHPC_VERSION/comm_libs/mpi/man

export LD_LIBRARY_PATH=$NVCOMPILERS/$NVHPC_VERSION/compilers/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NVCOMPILERS/$NVHPC_VERSION/comm_libs/mpi/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NVCOMPILERS/$NVHPC_VERSION/compilers/extras/qd/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$NVCOMPILERS/$NVHPC_VERSION/cuda/12.4/targets/x86_64-linux/lib:$LD_LIBRARY_PATH

mpirun -np 1 --map-by numa:PE=8 --bind-to core --report-bindings \
  -x OMP_NUM_THREADS=8 \
  -x OMP_PLACES=cores \
  -x OMP_PROC_BIND=close \
  --mca coll_hcoll_enable 0 \
  /home15/liangqm/ChihayaAnon-main/ShengBTE
```
