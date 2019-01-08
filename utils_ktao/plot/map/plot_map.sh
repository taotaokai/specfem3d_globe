#!/bin/bash

# plot s362ani topo map

#====== plot merged topo
gmt gmtset FONT_TITLE 16
gmt gmtset MAP_TITLE_OFFSET 0p
gmt gmtset MAP_FRAME_TYPE FANCY 
gmt gmtset PS_PAGE_ORIENTATION portrait
gmt gmtset MAP_GRID_PEN_PRIMARY 0.1p

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

out_fig=event_map.ps
psmeca_input="event_used.psmeca"

# basemap
gmt psbasemap -J$J -R$R -BWeSn+t"event map" -Ba5f1 -Xc0.0 -Yc0.0 -K -V > $out_fig

# topography
gmt grdimage topo.grd -R -J -V -O -K -Cgeo >> $out_fig

# coastal lines
gmt pscoast -J -R -W1.5p -N1/0.5p -Df -A1000 -O -K -V >> $out_fig

# geological blocks
gmt psxy zhangpz_block.txt -J -R -W1.0p,white -O -K -V >> $out_fig

# mesh boundaries
#awk '{print $2,$1}' mesh_corner_lla.txt |\
#  gmt psxy -J -R -W1.5p,blue -G- -O -K -V >> $out_fig

# stations
awk '{print $4,$3}' STATIONS | gmt psxy -J -R -St0.3c -Gblue -W0.5p,black -O -K -V >> $out_fig

# focal mechanisms
awk -F'|' '{printf "./make_psmeca_input.sh %s GCMT_1976-NOW.ndk\n",$9}' event_used.txt | bash > event_used.psmeca
cat $psmeca_input | gmt psmeca -J -R -Sm0.4 -Gred -L0.5p,black -W1.0p,white -O -V >> $out_fig
#awk '$3<=30' $psmeca_input | gmt psmeca -J -R -Sm0.4 -Gred -O -K -V >> $out_fig
#awk '$3>30&&$3<=120' $psmeca_input | gmt psmeca -J -R -Sm0.4 -Ggreen -O -K -V >> $out_fig
#awk '$3>120' $psmeca_input | gmt psmeca -J -R -Sm0.4 -Gblue -O -V >> $out_fig
#gmt makecpt -Cjet -T0/200/1  > my.cpt
#gmt psmeca $psmeca_input  -J -R -Sm0.4 -Zmy.cpt -O -K -V >> $out_fig
# colorbar
#gmt psscale -DjBM+w10c/0.5c+h -Y-8c -Bx50+l"event depth (km)" -Cmy.cpt -R -J -O -V >> $out_fig

ps2pdf $out_fig
evince ${out_fig//ps/pdf}
