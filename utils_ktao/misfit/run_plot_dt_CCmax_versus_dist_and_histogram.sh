#!/bin/bash

event_list=${1?[arg]need event list}
 
#for event_id in $(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
#do
#  echo ====== $event_id
#  for window_id in surface_10-25_Z surface_20-40_Z surface_35-60_Z surface_50-100_Z
#  do
#    echo ====== $window_id
#    ./plot_hist.py ${event_id}.txt ${window_id} ${event_id}_${window_id}.pdf
#  done
#done

#event_id=C201809080231A
#window_id=surface_50-100_Z
#window_id=surface_35-60_Z

#$(awk -F"|" 'NF&&$1!~/#/{printf "%s.%s.%s.%s\n", $1,$2,$3,$4}' $event_list)


for model in Shen2016 FWEA18 Xin2018 
do

  cat /dev/null > combined_${model}.txt
  for event_id in $(awk -F"|" 'NF&&$1!~/#/{printf "%s.%s.%s.%s\n", $1,$2,$3,$4}' $event_list) #$(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
  do
    cat ${model}/misfit_hist_noiseCC/${event_id}.txt >> combined_${model}.txt
  done

  #for window_id in p,P_Z s,S_Z s,S_T surface_10-25_Z surface_20-40_Z surface_35-60_Z # surface_50-100_Z
  for window_id in surface_10-25_Z surface_20-40_Z surface_35-60_Z # surface_50-100_Z
  do
 
    echo ====== combined $model $window_id
    ./plot_dt_CCmax_versus_dist_and_histogram.py combined_${model}.txt ${window_id} combined_${model}_${window_id}.pdf "${model} ${window_id}"

    for event_id in $(awk -F"|" 'NF&&$1!~/#/{printf "%s.%s.%s.%s\n", $1,$2,$3,$4}' $event_list) #$(awk -F"|" 'NF&&$1!~/#/{print $9}' $event_list)
    do
      echo ====== $model $window_id $event_id 
      ./plot_dt_CCmax_versus_dist_and_histogram.py ${model}/misfit_hist_noiseCC/${event_id}.txt ${window_id} ${model}_${event_id}_${window_id}.pdf  "${model} ${event_id}"
    done

  done

done