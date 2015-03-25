# -------------------------------------------------------------------------------- 
#    plotSchwGeods.gnu
#  
#  Copyright (c) 2009  Thomas Mueller, Frank Grave
#
#   This file is part of the m4d-library.
# 
# -------------------------------------------------------------------------------- 
#     Run e.g.:  ./m4dTestGeodesic 
# -------------------------------------------------------------------------------- 

#set term fig
#set output "schw_geods.fig"

set term png
set output "schw_geods.png"

set xrange[-10:10]
set yrange[-10:10]
set size ratio 1
unset grid
unset key

set parametric
plot [t=0:2*pi] "points_0.dat" u 3:4 w l lt 1, "points_1.dat" u 3:4 w l lt 1, "points_2.dat" u 3:4 w l lt 1, "points_3.dat" u 3:4 w l lt 1, "points_4.dat" u 3:4 w l lt 1, "points_5.dat" u 3:4 w l lt 1, "points_6.dat" u 3:4 w l lt 1, "points_7.dat" u 3:4 w l lt 1, "points_8.dat" u 3:4 w l lt 1, "points_9.dat" u 3:4 w l lt 1, "points_10.dat" u 3:4 w l lt 1, "points_11.dat" u 3:4 w l lt 1, "points_12.dat" u 3:4 w l lt 1, "points_13.dat" u 3:4 w l lt 1,"points_14.dat" u 3:4 w l lt 1, 2*cos(t),2*sin(t) lt 3

