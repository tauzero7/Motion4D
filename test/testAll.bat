@ rem ------------------------
@ rem  Run all test programs
@ rem ------------------------
@ SET CONFIG=Release
@ IF "%1"=="Debug" SET CONFIG=Debug
@ echo Run %CONFIG% test programs

@ SET DAR=
@ IF %CONFIG%==Debug  SET DAR=_d

@ cd testDatabase
.\m4dTestDatabase%DAR%.exe
@ cd ..

@ cd calcGeodesic
.\m4dCalcGeodesic%DAR%.exe kerr.ini  > points_kerr.dat
@ cd ..

@ cd calcParallel
.\m4dCalcParallel%DAR%.exe kerr_timelike.ini > points_par.dat
@ cd ..

@ cd testCircular
.\m4dTestCircular%DAR%.exe 6  > points_circ.dat 
@ cd ..

@ cd testFWworldline
.\m4dTestFWworldline%DAR%.exe > points_fw.dat
@ cd ..

@ cd testGeodesic
.\m4dTestGeodesic%DAR%.exe
@ cd ..

@ cd testJacobi
.\m4dTestJacobi%DAR%.exe kerr_lens.ini > points_lens.dat
@ cd ..

