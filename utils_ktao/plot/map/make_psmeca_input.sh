#!/bin/bash

# create gmt psmeca inputs from GCMT ndk file for the given event gcmt-ID

# -S
# (m) Seismic moment tensor (Harvard CMT, with zero trace):
#      X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, newX, newY, event_title

#====== command line args
evid=${1:?[arg] need evid}
ndk_file=${2:?[arg] need GCMT_1976-NOW.ndk}
#outfile=${3:?[arg] need out filename}

#====== 
tmp_ndk=$(mktemp)
grep -A 3 -B 1 "^$evid" $ndk_file > $tmp_ndk

#sed '1!d;s/[\/:]/ /g;q' $tmp_ndk  > $outfile

# X, Y, depth,
sed '3!d;q' $tmp_ndk | awk '{printf "%8.4f %8.4f %8.4f ",$6,$4,$8}'

# mrr, mtt, mff, mrt, mrf, mtf, exp,
sed '4!d;q' $tmp_ndk | awk '{printf "%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %02d ",$2,$4,$6,$8,$10,$12,$1}'

## newX, newY
#sed '2!d;q' $tmp_ndk | awk '{printf "0.0 0.0\n"}'

# newX, newY, event_title
sed '2!d;q' $tmp_ndk | awk '{printf "0.0 0.0 %s\n",$1}'

rm $tmp_ndk
