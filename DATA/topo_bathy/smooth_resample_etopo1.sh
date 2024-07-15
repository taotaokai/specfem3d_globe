#!/bin/bash

# ktao: etopo1 is smoothed and resampled

input=ETOPO1_Ice_g_gmt4.grd
output=ETOPO1_Ice_g_smooth_b15km_I4m.grd

gmt grdfilter $input -G${output} -D2 -Fb15 -I4m
