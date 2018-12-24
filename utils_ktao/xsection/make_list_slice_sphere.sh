#!/bin/bash

# make a list of spherical slices for xsem_slice_sphere

out_file=${1:?[arg]need output list (e.g. slice_sphere.list)}

#---- geographic range and mesh size
lat0=21
lat1=35
nlat=29
lon0=94
lon1=110
nlon=33

#---- depths(km): 5, 10, 20:20:1000 
#echo 3 > depth.list
#seq 5 5 900 > depth.list

#---- make list
echo "# spherical slices:" > $out_file

#cat depth.list |\
#while read depth
for depth in 3 5 15 20 40 60 80 100 120
do
  printf "%6.2f %6.2f %4d %7.2f %7.2f %4d %6.1f depth_%03d\n" \
    $lat0 $lat1 $nlat $lon0 $lon1 $nlon $depth $depth
done >> $out_file
