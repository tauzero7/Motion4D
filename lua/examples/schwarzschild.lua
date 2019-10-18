--[[
  Test file for calculating geodesic in the Schwarzschild spacetime
]]

--PrintObject();
--PrintMetricDB();

SetMetric("Schwarzschild",{mass = 1.0})
--SetMetricParams({mass = 2.0})
PrintMetric()

PrintSolverDB()

SetGeodSolver("GSL_RK_Fehlberg")
SetGeodParams({stepctrl = "yes",
               stepsize = 1e-2,
               epsabs = 1e-8,
               constr = 1e-6})
PrintSolver()


SetLocalTetrad({
    pos = {0.0, 10.0, math.pi*0.5, 0.0},
    type = 'default',
    e0 = {1.0, 0.0, 0.0, 0.0},
    e1 = {0.0, 0.0, 0.0, 1.0},
    e4 = {0.0, 0.0,-1.0, 0.0},
    e2 = {0.0,-1.0, 0.0, 0.0}
})
PrintTetrad()

phi = math.rad(30.0)
numPoints = CalculateGeodesic({
    dir = {math.sin(phi), math.cos(phi), 0.0},
    type = "lightlike",
    timedir = "backward",
    maxpoints = 30
})

--SaveGeodesic("schw_geod.dat")

fptr = io.open("schw_geod.dat","w")
if fptr == nil then
    print("Cannot open file!\n")
end

for i=0, numPoints-1, 1 do    
    pos = GetGeodPoint(i)    
    --PrintVec(pos,"%14.6f ");    
    newpos = CoordTrans(pos,"spherical","cartesian")
    PrintVec(newpos,"%14.6f ");
    --fptr.write(fptr,string.format("%14.6f %14.6f %14.6f %14.6f\n",newpos[1],newpos[2],newpos[3],newpos[4]))
end

io.close(fptr)


mat = CalcTidalMatrix({pos = {0.0, 3.0, math.pi*0.5, 0.0}})
PrintMat(mat, "%12.8f ")
