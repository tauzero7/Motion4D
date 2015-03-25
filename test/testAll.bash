#!/bin/bash

DR=""
if [ "$1" == "debug" ]; then
	DR="_debug"
fi

if [ -d testDatabase ]; then
    cd testDatabase
    ./m4dTestDatabase$DR
    cd ..
else
	./m4dTestDatabase$DR
fi



if [ -d calcGeodesic ]; then
	cd calcGeodesic
    ./m4dCalcGeodesic$DR  kerr.ini  > points_kerr.dat
    gnuplot plotGeodesic.gnu
	cd ..
else
	./m4dCalcGeodesic$DR  kerr.ini  > points_kerr.dat
	gnuplot plotGeodesic.gnu
fi



if [ -d calcParallel ]; then
	cd calcParallel
    ./m4dCalcParallel$DR kerr_timelike.ini > points_par.dat
	cd ..
else
    ./m4dCalcParallel$DR kerr_timelike.ini > points_par.dat
fi



if [ -d testCircular ]; then
	cd testCircular
    ./m4dTestCircular$DR 6  >  points_circ.dat 
    cd ..
else
    ./m4dTestCircular$DR 6  >  points_circ.dat 
fi



if [ -d testFWworldline ]; then
	cd testFWworldline
    ./m4dTestFWworldline$DR > points_fw.dat
    cd ..
else
    ./m4dTestFWworldline$DR > points_fw.dat
fi


if [ -d testGeodesic ]; then
	cd testGeodesic
    ./m4dTestGeodesic$DR
	gnuplot plotSchwGeods.gnu
	cd ..
else
    ./m4dTestGeodesic$DR
	gnuplot plotSchwGeods.gnu
fi



if [ -d testJacobi ]; then
	cd testJacobi
    ./m4dTestJacobi$DR kerr_lens.ini > points_lens.dat
	gnuplot plotJacobi.gnu
	cd ..
else
	./m4dTestJacobi$DR kerr_lens.ini > points_lens.dat
	gnuplot plotJacobi.gnu
fi



