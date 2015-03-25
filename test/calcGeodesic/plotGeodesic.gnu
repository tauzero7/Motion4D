# -------------------------------------------------------------------------------- 
#    plotGeodesic.gnu
#  
#  Copyright (c) 2009  Thomas Mueller, Frank Grave
#
#   This file is part of the m4d-library.
# 
# -------------------------------------------------------------------------------- 
#     Run e.g.:  ./m4dCalcGeodesic kerr.ini  > points_kerr.dat
# -------------------------------------------------------------------------------- 

#set term fig
#set output "kerr_geod.fig"

set term png
set output "kerr_geod.png"

set xrange[-3:3]
set yrange[-3:3]
set size ratio 1
set grid
unset key

set ticslevel 0
set view 70,34
set size 1.2,0.94

splot "points_kerr.dat" u 3:4:5 w l

