@ rem ------------------------
@ rem  Run all test programs
@ rem ------------------------
@ SET CONFIG=Debug
@ IF "%1"=="Release" SET CONFIG=Release
@ echo Run %CONFIG% test programs

@ SET DAR=
@ IF %CONFIG%==Debug  SET DAR=d

@ cd testDatabase
.\%CONFIG%\m4dTestDatabase%DAR%.exe
@ cd ..

@ cd calcGeodesic
.\%CONFIG%\m4dCalcGeodesic%DAR%.exe kerr.ini  > points_kerr.dat
@ cd ..

@ cd calcParallel
.\%CONFIG%\m4dCalcParallel%DAR%.exe kerr_timelike.ini > points_par.dat
@ cd ..

@ cd testCircular
.\%CONFIG%\m4dTestCircular%DAR%.exe 6  > points_circ.dat 
@ cd ..

@ cd testFWworldline
.\%CONFIG%\m4dTestFWworldline%DAR%.exe > points_fw.dat
@ cd ..

@ cd testGeodesic
.\%CONFIG%\m4dTestGeodesic%DAR%.exe
@ cd ..

@ cd testJacobi
.\%CONFIG%\m4dTestJacobi%DAR%.exe kerr_lens.ini > points_lens.dat
@ cd ..

