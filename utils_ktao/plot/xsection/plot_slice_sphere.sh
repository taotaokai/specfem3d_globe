#!/bin/bash

# plot xsections for model relative difference
sem_utils_dir=${1:?[arg]need sem_utils_dir for plot_slice_sphere.py}
nc_dir=${2:?[arg]need model_dir (*.nc)}
slice_list=${3:?[arg]need slice_list}
out_dir=${4:?[arg]need out_dir}
model_names=${5?[arg]need model_names, e.g. vsv,vpv}
python_exec=${6?[arg]need python exec}

awk 'NF&&$1!~/#/' $slice_list |\
while read lat0 lat1 nlat lon0 lon1 nlon depth nc_tag 
do

  echo
  echo "#====== $nc_tag: $lat0 $lat1 $nlat $lon0 $lon1 $nlon $depth"
  echo

  $python_exec $sem_utils_dir/plot_slice_sphere.py  $nc_dir/${nc_tag}.nc $model_names  $out_dir/${nc_tag}.pdf

done