#!/bin/bash

# interpolate SEM gll models to create spherical slices

#====== command line args
sem_utils_dir=${1:?[arg]need sem_utils_dir dir to find sem_slice_gcircle}
mesh_dir=${2:?[arg]need mesh_dir for proc*_reg1_solver_data.bin}
model_dir=${3:?[arg]need model_dir for proc*_reg1_<model_name>.bin}
slice_list=${4:?[arg]need slice_list}
model_names=${5:?[arg]need model names, e.g. vsv,vsh,vpv,vph,rho,eta}
model_units=${6:?[arg]need model units, e.g. km.s-1,g.cm-3,1}
nc_dir=${7:?[arg]need output directory for .nc files}
job_file=${8:?[arg]need output job_file}
mpi_exec=${9:?[arg]need mpi_exec, e.g. ibrun or mpirun}

python_exec=mypython

# check input parameters
sem_utils_dir=$(readlink -f $sem_utils_dir)
if [ ! -d "$sem_utils_dir" ]
then
  echo "[ERROR] sem_utils_dir not found: " $sem_utils
  exit -1
fi

mesh_dir=$(readlink -f $mesh_dir)
if [ ! -d "$mesh_dir" ]
then
  echo "[ERROR] mesh_dir not found: " $mesh_dir
  exit -1
fi

model_dir=$(readlink -f $model_dir)
if [ ! -d "$model_dir" ]
then
  echo "[ERROR] model_dir not found: " $model_dir
  exit -1
fi

slice_list=$(readlink -f $slice_list)
if [ ! -f "$slice_list" ]
then
  echo "[ERROR] slice_list not found: " $slice_list
  exit -1
fi

# get some mesh info
nproc=$(ls $mesh_dir/proc*_reg1_solver_data.bin | wc -l)

#====== create job script
cat <<EOF > $job_file
#!/bin/bash
#SBATCH -J slice_sphere
#SBATCH -o ${job_file}.o%j
#SBATCH -N 1
##SBATCH -n $nproc
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -t 09:00:00
##SBATCH --mail-user=kai.tao@utexas.edu
##SBATCH --mail-type=begin
##SBATCH --mail-type=end

echo
echo "Start: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo

mkdir $nc_dir

EOF

grep -v "^#" $slice_list |\
while read lat0 lat1 nlat lon0 lon1 nlon depth fname
do

cat<<EOF >> $job_file
echo
echo \$(date)
echo "#====== $fname: $depth "
echo
#${mpi_exec}
${python_exec} \
    $sem_utils_dir/sem_slice_sphere.py \
    $nproc \
    $mesh_dir \
    $model_dir \
    $model_names \
    $model_units \
    $lat0 $lat1 $nlat \
    $lon0 $lon1 $nlon \
    $depth \
    $nc_dir/${fname}.nc

EOF

done

cat<<EOF >> $job_file
echo
echo "End: JOB_ID=\${SLURM_JOB_ID} [\$(date)]"
echo
EOF
