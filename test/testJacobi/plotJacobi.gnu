# -------------------------------------------------------------------------------- 
#    plotJacobi.gnu
#  
#  Copyright (c) 2009,2010  Thomas Mueller, Frank Grave
#
#   This file is part of the m4d-library.
# 
# -------------------------------------------------------------------------------- 
#     Run:  ./m4dTestJacobi  kerr_lens.ini -> points_lens.dat
# -------------------------------------------------------------------------------- 

#set term fig
#set output "kerr_lens.fig"

set term png
set output "kerr_lens.png"

set xrange[0:1800]
set yrange[0.0001:1000]
set logscale y
unset grid
unset key

plot "points_lens.dat" u 1:9 w l lt 1

