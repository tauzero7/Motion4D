--[[
  Calculate motion of a particle on the last stable circular orbit around
  a Schwarzschild black hole.
--]]

SetMetric("Schwarzschild", {mass = 1.0})
PrintMetric()

SetGeodSolver("GSL_RK_Fehlberg")
SetGeodParams({stepctrl = "yes",
               stepsize = 1e-2,
               epsabs = 1e-12,
               constr = 1e-6})
PrintSolver()

SetLocalTetrad({
    pos = {0.0, 6.0, math.pi*0.5, 0.0},
    type = 'default',
    e0 = {1.0, 0.0, 0.0, 0.0},
    e1 = {0.0, 0.0, 0.0, 1.0},
    e4 = {0.0, 0.0,-1.0, 0.0},
    e2 = {0.0,-1.0, 0.0, 0.0}
})
PrintTetrad()

numPoints = CalculateParTransport({
    dir = {0.0, 0.0, 1.0},
    beta = 0.5,
    timedir = 'forward',
    maxpoints = 1000
})

print("Write 'circOrbit.dat' ... ")
fptr = io.open("circOrbit.dat","w")
if fptr == nil then
    print("Cannot open file!\n")
end

for i=0, numPoints-1, 1 do    
    pos = GetGeodPoint(i)        
    newpos = CoordTrans(pos,"spherical","cartesian")
    --coord_e1 = GetTrajE1(i)
    --local_e1 = CoordToLocal({pos = pos, vec = coord_e1})
    --PrintVec(local_e1)
    --phi = -local_e1[1]^2 + local_e1[2]^2 + local_e1[3]^2 + local_e1[4]^2
    --print(phi)
        
    fptr.write(fptr,string.format("%14.6f %14.6f %14.6f %14.6f\n",
        newpos[1],newpos[2],newpos[3],newpos[4]))
end

io.close(fptr)
