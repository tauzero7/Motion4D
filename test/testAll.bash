#!/bin/bash

PYTHON=ipython

DR=""
if [ "$1" == "debug" ]; then
	DR="_d"
fi


if [ -d testDatabase ]; then
    cd testDatabase
    ./m4dTestDatabase$DR
    cd ..
fi


if [ -d calcGeodesic ]; then
    cd calcGeodesic
    ./m4dCalcGeodesic$DR  kerr.ini  > points_kerr.dat
    #gnuplot plotGeodesic.gnu
    $PYTHON plotGeodesic.py
    cd ..
fi


if [ -d calcParallel ]; then
    cd calcParallel
    ./m4dCalcParallel$DR kerr_timelike.ini > points_par.dat
    cd ..
fi


if [ -d testCircular ]; then
    cd testCircular
    ./m4dTestCircular$DR 6  >  points_circ.dat 
    cd ..
fi


if [ -d testFWworldline ]; then
    cd testFWworldline
    ./m4dTestFWworldline$DR > points_fw.dat
    cd ..
fi


if [ -d testGeodesic ]; then
    cd testGeodesic
    ./m4dTestGeodesic$DR
    #gnuplot plotSchwGeods.gnu
    $PYTHON plotSchwGeods.py
    cd ..
fi


if [ -d testJacobi ]; then
    cd testJacobi
    ./m4dTestJacobi$DR kerr_lens.ini > points_lens.dat
    #gnuplot plotJacobi.gnu
    $PYTHON plotJacobi.py
    cd ..

fi

