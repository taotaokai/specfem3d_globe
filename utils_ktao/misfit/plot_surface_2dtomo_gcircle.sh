#!/bin/bash

#====== plot merged topo
gmt set \
  FONT 22p,Times-Roman \
  FONT_TITLE 18 \
  MAP_TITLE_OFFSET 0p \
  MAP_FRAME_TYPE FANCY \
  PS_PAGE_ORIENTATION portrait \
  MAP_GRID_PEN_PRIMARY 0.1p

nc_file=${1:?[arg]need nc_file}
out_fig=${2:?[arg]need out_fig}
plot_title=${3:?[arg]need plot_title}
#out_fig=${nc_file/.nc/.ps}

#nc_file=combined_Shen2016_surface_10-25_Z.nc
#out_fig=event_map.ps

lon0=93
lon1=110
lat0=21
lat1=35
map_size=15.0c

lon_mid=$(echo "($lon0 + $lon1)/2" | bc -l)
lat_mid=$(echo "($lat0 + $lat1)/2" | bc -l)

R=${lon0}/${lon1}/${lat0}/${lat1}
#J=T${lon_mid}/${lat_mid}/${map_size}
J=M${lon_mid}/${map_size}

# basemap
gmt psbasemap -J$J -R$R -BWeSn+t"$plot_title" -Ba5f1 -Xc-1c -Yc0.0 -K -V > $out_fig

# topography
#gmt grdimage topo.grd -R -J -V -O -K -Cgeo >> $out_fig
gmt makecpt -T-0.02/0.02/0.001 -N -Z -Cjet > ds.cpt

gmt grdmath ${nc_file}?path_density NORM 0.5 SUB = mask.nc
gmt grdimage ${nc_file}?ds -Imask.nc -R -J -Cds.cpt -O -K -V >> $out_fig
gmt grdcontour ${nc_file}?path_density -R -J -C+0.1 -W2p,black,solid -O -K -V >> $out_fig
#gmt psscale -DJRM+w8c/0.3c+o1c/0c+e -By+l"  (s/km)" -Cds.cpt -Baf -R -J -O -K -V >> $out_fig

# coastal lines
gmt pscoast -J -R -W1.5p -N1/0.5p -Df -A1000 -O -K -V >> $out_fig

# geological blocks
gmt psxy zhangpz_block.txt -J -R -W1.0p,white -O -K -V >> $out_fig

## stations for noiseCC
#awk '{print $4,$3}' STATIONS | gmt psxy -J -R -Ss0.2c -Gblack -O -K -V >> $out_fig
## virtual sources for noiseCC
#awk -F"|" '{print $6,$5}' event_noiseCC.txt | gmt psxy -J -R -Sa0.5c -Gred -W0.5p,black -O -K -V >> $out_fig

# stations for regEQ
awk '{print $4,$3}' STATIONS_merge | gmt psxy -J -R -Ss0.2c -Gblack -O -K -V >> $out_fig
# sources for regEQ 
cat event_regEQ.psmeca | gmt psmeca -J -R -Sm0.4 -Gblack -L0.5p,black -W1.0p,white -O -K -V >> $out_fig

# insert histogram
gmt psbasemap -R -J -DjBL+w5.0c/2.7c+o0.0c/0.0c+stmp -F+gwhite -O -K -V >> $out_fig
#read x0 y0 w h < tmp
gmt set MAP_TICK_LENGTH -0.1c FONT_ANNOT 10p MAP_FRAME_PEN 0.5p

ncdump -v dt_res -f c ${nc_file} | sed -e '1,/data:/d' -e '$d' | sed "s/,.*//;s/;.*//;s/.*=//;s/[ ]*//" |\
  gmt pshistogram -JX4c/2c -R-5/5/0/500 -BWSne -Bx5.0 -By300 -W0.5p -Gred -X0.8c -Y0.4c -F -O -K -V >> $out_fig

ncdump -v dt -f c ${nc_file} | sed -e '1,/data:/d' -e '$d' | sed "s/,.*//;s/;.*//;s/.*=//;s/[ ]*//" |\
  gmt pshistogram  -J -R -W0.5p -L1p,black -F -O -V >> $out_fig

ps2pdf $out_fig
#evince ${out_fig//ps/pdf}
