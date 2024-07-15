#!/bin/bash

for misfit_file in *.LHZ.txt
do

  echo "#------ $misfit_file"

  pattern=$(echo $misfit_file | awk -F. '{printf "%s[ ]*%s.", $1,$2}')
  evlo=$(grep "$pattern" STATIONS | awk '{print $4}')
  evla=$(grep "$pattern" STATIONS | awk '{print $3}')

  while read -u 10 net sta stla stlo junk
  do
    
    pattern="${net}.${sta}"
    grep "^$pattern" $misfit_file | sed "s/$/ |  $stlo $stla $evlo $evla/"

  done 10< STATIONS 

done
