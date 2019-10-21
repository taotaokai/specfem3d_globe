#!/bin/bash

# create SEM build directory and compile the program

sem_source_dir=${1:?[arg] need sem_source_dir}
sem_setup_dir=${2:?[arg] need setup dir(for constants.h.in)}
sem_data_dir=${3:?[arg] need DATA dir(for DATA/Par_file,STATIONS,CMTSOLUTION)}
sem_build_dir=${4:?[arg] need build dir(e.g. specfem3d_globe)}

#mpif90=${5:?[arg]need mpif90 executable}
#mpicc=${6:?[arg]need mpicc executable}

sem_source_dir=$(readlink -f $sem_source_dir)
sem_setup_dir=$(readlink -f $sem_setup_dir)
sem_data_dir=$(readlink -f $sem_data_dir)
sem_build_dir=$(readlink -f $sem_build_dir)

if [ x${sem_source_dir} = x${sem_build_dir} ]
then
  echo "[ERROR] sem_source_dir cannot be the same as sem_build_dir!"
  exit -1
fi
if [ ! -d "$sem_source_dir" ]
then
  echo "[ERROR] sem_source_dir does not exits!"
  exit 1
fi
if [ ! -d "$sem_setup_dir" ]
then
  echo "[ERROR] sem_setup_dir does not exits!"
  exit 1
fi
if [ ! -d "$sem_data_dir" ]
then
  echo "[ERROR] sem_data_dir does not exits!"
  exit 1
fi
if [ ! -f "$sem_data_dir/CMTSOLUTION" ]
then
  echo "[ERROR] sem_data_dir/DATA/CMTSOLUTION does not exits!"
  exit 1
fi

#====== prepare sem_build_dir
if [ -d "$sem_build_dir" ]
then
  echo "[WARN] sem_build_dir exits, delete!"
  rm -rf $sem_build_dir
fi
mkdir $sem_build_dir

cd $sem_build_dir
mkdir DATA DATABASES_MPI OUTPUT_FILES

# link src/, config*
ln -s ${sem_source_dir}/src ./
ln -s ${sem_source_dir}/config.guess ./
ln -s ${sem_source_dir}/config.sub ./
ln -s ${sem_source_dir}/configure ./
ln -s ${sem_source_dir}/configure.ac ./
ln -s ${sem_source_dir}/flags.guess ./
ln -s ${sem_source_dir}/install-sh ./
ln -s ${sem_source_dir}/Makefile.in ./
# link DATA/*
ln -s ${sem_source_dir}/DATA/* ${sem_build_dir}/DATA/

# copy setup directory
cp -r ${sem_source_dir}/setup ${sem_build_dir}/
chmod u+w -R setup

#====== use sem_config_dir
cd $sem_build_dir/DATA
#rm -f Par_file CMTSOLUTION
cp $sem_data_dir/Par_file .
cp $sem_data_dir/CMTSOLUTION .
#cp $sem_data_dir/FORCESOLUTION .
cp $sem_data_dir/STATIONS .

cd $sem_build_dir/setup
rm -f constants.h.in
cp $sem_setup_dir/constants.h.in .

#====== build SEM
cd $sem_build_dir

#./configure FC=gfortran MPIFC=mpif90 FCFLAGS="-O3 -lpthread" CC=mpicc CFLAGS="-O3"
#./configure FC=gfortran MPIFC=mpif90 FCFLAGS="-march=core-avx2 -O3 -lpthread" CC=mpicc CFLAGS="-O3 -march=core-avx2"
#./configure FC=gfortran CC=gcc MPIFC=${mpif90} FCFLAGS="-O3 -lpthread" CFLAGS="-O3"
./configure FC=gfortran CC=gcc MPIFC=mpif90 MPILIBS="-lpthread"

make clean
#make xdecompose_mesh xmeshfem3D xgenerate_databases xspecfem3D xcombine_vol_data_vtk \
#  > >(tee make.out) 2> >(tee make.err >&2)
make all  > >(tee make.out) 2> >(tee make.err >&2)

#END