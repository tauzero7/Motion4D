/* -------------------------------------------------------------------------------
    m4dlua_solver.cpp

  Copyright (c) 2015  Thomas Mueller

   This file is part of the m4d-library.

   The m4d-library is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The m4d-library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with the m4d-library.  If not, see <http://www.gnu.org/licenses/>.
  ------------------------------------------------------------------------------- */

#include "lua/m4dlua_solver.h"
#include "math/TransCoordinates.h"

/**
 * @brief Register functions for lua.
 * @param L  Pointer to lua state.
 */
void lua_reg_solver(lua_State *L) {
    lua_register(L, "PrintSolverDB", printSolverDB);
    lua_register(L, "PrintGeodTypes", printGeodTypes);

    lua_register(L, "SetGeodSolver", setGeodSolver);
    lua_register(L, "SetGeodParams", setGeodParams);
    lua_register(L, "PrintSolver", printSolver);

    lua_register(L, "SetLocalTetrad", setLocalTetrad);
    lua_register(L, "PrintTetrad", printTetrad);

    lua_register(L, "CalculateGeodesic", calculateGeodesic);

    lua_register(L, "InterpPosByCoord", interpolatePosByCoord);
    lua_register(L, "InterpPosByLambda", interpolatePosByLambda);

    lua_register(L, "GetCurrentDir", getCurrentDir);
    lua_register(L, "GetCurrentPos", getCurrentPos);
    lua_register(L, "GetCurrentLambda", getCurrentLambda);

    lua_register(L, "SaveGeodesic", saveGeodesic);
}


/**
 * @brief Print database of all solvers to stdout.
 * @return 0 : no arguments will be put onto lua stack
 */
int printSolverDB(lua_State* ) {
    for(unsigned int i=0; i<m4d::NUM_GEOD_SOLVERS; i++) {
        fprintf(stdout,"%s\n",m4d::stl_solver_nicknames[i]);
    }
    return 0;
}


/**
 * @brief Print all geodesic types to stdout.
 * @return 0 : no arguments will be put onto lua stack
 */
int printGeodTypes(lua_State* ) {
    for(unsigned int i=0; i<m4d::NUM_ENUM_GEODESIC_TYPE; i++)   {
        fprintf(stdout,"%2d .. %s\n",i,m4d::stl_geodesic_type[i]);
    }
    return 0;
}


/**
 * @brief Set geodesic solver.
 * @param L Pointer to lua state.
 * @return 0 : no arguments will be put onto lua stack
 */
int setGeodSolver(lua_State* L) {
    int numParams = lua_gettop(L);
    if (numParams==0) {
        mlua_error(L,"SetGeodesicSolver: missing solver name!\n");
        return 0;
    }

    if (mObject.currMetric==NULL) {
        mlua_error(L,"SetGeodesicSolver: missing metric!\n");
        return 0;
    }

    if (!lua_isnil(L,-numParams)) {
        if (lua_isstring(L,-numParams)) {
            std::string solverName = std::string(lua_tostring(L,-numParams));
            //std::cerr << solverName << std::endl;
            mObject.geodSolver = mObject.solverDB->getIntegrator(mObject.currMetric,solverName);
            mObject.geodSolverType = mObject.solverDB->getIntegratorNr(solverName);
        }
    }
    lua_pop(L,lua_gettop(L));
    return 0;
}


/**
 * @brief Set geodesic parameters.
 *
 *   The parameters have to be given in table format with key-value pairs.
 * @param L Pointer to lua state.
 * @return 0 : no arguments will be put onto lua stack
 */
int setGeodParams(lua_State* L) {
    if (mObject.geodSolver == NULL) {
        mlua_error(L,"SetGeodParams: you first have to set a geodesic solver!\n");
        return 0;
    }

    std::string sVal;
    double      pVal;

    if (!lua_isnil(L,-1)) {
        if (lua_istable(L,-1)) {
            if (getfield(L,"stepctrl",sVal)) {
                if (sVal.compare("yes")==0) {
                    mObject.stepsizeControlled = true;
                    mObject.geodSolver->setStepSizeControlled(true);
                } else {
                    mObject.stepsizeControlled = false;
                    mObject.geodSolver->setStepSizeControlled(false);
                }
            }

            if (getfield(L,"stepsize",pVal)) {
                mObject.stepsize = pVal;
                mObject.geodSolver->setAffineParamStep(pVal);
            }

            if (getfield(L,"epsabs",pVal)) {
                double eps_a,eps_r;
                mObject.epsAbs = pVal;
                mObject.geodSolver->getEpsilons(eps_a,eps_r);
                mObject.geodSolver->setEpsilons(pVal,eps_r);
            }

            if (getfield(L,"epsrel",pVal)) {
                double eps_a,eps_r;
                mObject.epsRel = pVal;
                mObject.geodSolver->getEpsilons(eps_a,eps_r);
                mObject.geodSolver->setEpsilons(eps_a,pVal);
            }

            if (getfield(L,"constr",pVal)) {
                mObject.epsConstr = pVal;
                mObject.geodSolver->setConstrEps(pVal);
            }

            if (getfield(L,"type",sVal)) {
                if (sVal.compare("lightlike")==0 || sVal.compare("null")==0) {
                    mObject.type = m4d::enum_geodesic_lightlike;
                }
                else if (sVal.compare("spacelike")==0) {
                    mObject.type = m4d::enum_geodesic_spacelike;
                }
                else if (sVal.compare("timelike")==0) {
                    mObject.type = m4d::enum_geodesic_timelike;
                }
                else if (sVal.compare("lightlikeSachs")==0) {
                    mObject.type = m4d::enum_geodesic_lightlike_sachs;
                }
                mObject.geodSolver->setGeodesicType(mObject.type);
            }

            if (getfield(L,"max_stepsize",pVal)) {
                mObject.max_stepsize = pVal;
                mObject.geodSolver->setMaxAffineParamStep(pVal);
            }

            if (getfield(L,"min_stepsize",pVal)) {
                mObject.min_stepsize = pVal;
                mObject.geodSolver->setMinAffineParamStep(pVal);
            }

            std::vector<double> bbmin;
            if (getfield(L, "bbmin", bbmin) && bbmin.size()>3) {
                m4d::vec4 bboxMin, bboxMax;
                mObject.geodSolver->getBoundingBox(bboxMin,bboxMax);
                bboxMin = m4d::vec4(bbmin[0], bbmin[1], bbmin[2], bbmin[3]);
                mObject.geodSolver->setBoundingBox(bboxMin,bboxMax);
            }

            std::vector<double> bbmax;
            if (getfield(L, "bbmax", bbmax) && bbmax.size()>3) {
                m4d::vec4 bboxMin, bboxMax;
                mObject.geodSolver->getBoundingBox(bboxMin,bboxMax);
                bboxMax = m4d::vec4(bbmax[0], bbmax[1], bbmax[2], bbmax[3]);
                mObject.geodSolver->setBoundingBox(bboxMin,bboxMax);
            }
        }
    }
    return 0;
}


/**
 * @brief Print solver data to stdout.
 * @return 0 : no arguments will be put onto lua stack
 */
int printSolver(lua_State* ) {
    if (mObject.geodSolver!=NULL) {
        mObject.geodSolver->printF(stdout);
    }
    return 0;
}


/**
 * @brief Set the local tetrad vectors (base[0..3]) of m4dObject.
 * @param L Pointer to lua state.
 * @return 0 : no arguments will be put onto lua stack
 */
int setLocalTetrad(lua_State* L) {

    if (lua_isnil(L,-1) ||!lua_istable(L,-1)) {
        mlua_error(L,"LocalTetrad expects a table of values.\n");
        return 0;
    }

    double pVal;
    std::string sVal;

    std::vector<double> initPos;
    if (getfield(L,"pos",initPos) && initPos.size()==4) {
        mObject.startPos = m4d::vec4(initPos[0],initPos[1],initPos[2],initPos[3]);
    }

    if (getfield(L,"type",pVal)) {
        mObject.tetradType = (m4d::enum_nat_tetrad_type)((int)pVal);
    }
    if (getfield(L,"type",sVal)) {
        unsigned int i=0;
        while(i < m4d::NUM_ENUM_NAT_TETRAD_TYPES) {
            if (sVal.compare(m4d::stl_nat_tetrad_types[i])==0) {
                mObject.tetradType = static_cast<m4d::enum_nat_tetrad_type>(i);
                // test if metric allows this tetrad
                break;
            }
            i++;
        }
    }

    std::vector<double> e0;
    if (getfield(L,"e0",e0) && e0.size()==4) {
        mObject.base[0] = m4d::vec4(e0[0],e0[1],e0[2],e0[3]);
    }

    std::vector<double> e1;
    if (getfield(L,"e1",e1) && e1.size()==4) {
        mObject.base[1] = m4d::vec4(e1[0],e1[1],e1[2],e1[3]);
    }

    std::vector<double> e2;
    if (getfield(L,"e2",e2) && e2.size()==4) {
        mObject.base[2] = m4d::vec4(e2[0],e2[1],e2[2],e2[3]);
    }

    std::vector<double> e3;
    if (getfield(L,"e3",e3) && e3.size()==4) {
        mObject.base[3] = m4d::vec4(e3[0],e3[1],e3[2],e3[3]);
    }

    return 0;
}


/**
 * @brief Print tetrad data to stdout.
 * @param L Pointer to lua state.
 * @return 0 : no arguments will be put onto lua stack
 */
int printTetrad(lua_State* ) {
    fprintf(stdout,"\nTetrad:\n-------\n");
    fprintf(stdout,"\tpos  : ");mObject.startPos.printS(stdout);
    //fprintf(stderr,"\tldir : ");mObject.startDir.printS();
    fprintf(stdout,"\ttype : %s\n",m4d::stl_nat_tetrad_types[mObject.tetradType]);
    fprintf(stdout,"\te0   : ");mObject.base[0].printS(stdout);
    fprintf(stdout,"\te1   : ");mObject.base[1].printS(stdout);
    fprintf(stdout,"\te2   : ");mObject.base[2].printS(stdout);
    fprintf(stdout,"\te3   : ");mObject.base[3].printS(stdout);
    fprintf(stdout,"\n");
    return 0;
}


/**
 * @brief Calculate geodesic.
 *
 *   Before calculating a geodesic, a metric as well as a geodesic solver have to be defined.
 *   By means of a lua table, the geodesic type (type: lightlike, timelike, spacelike), the time
 *   direction (timedir: forward, backward), the initial velocity (vel), the maximum number of
 *   points (maxpoints), and the initial 3-direction (dir) should be specified.
 * @param L Pointer to lua state.
 * @return  Number of calculated points.
 */
int calculateGeodesic(lua_State* L) {
    if (mObject.geodSolver == NULL) {
        mlua_error(L,"CalculateGeodesic: no geodesic solver set!\n");
        return 0;
    }

    double      pVal;
    std::string sVal;

    if (!lua_isnil(L,-1)) {
        if (lua_istable(L,-1)) {
            std::vector<double> initDir;
            getfield(L,"dir",initDir);
            if (initDir.size()==3) {
                mObject.startDir = m4d::vec3(initDir[0],initDir[1],initDir[2]).getNormalized();
                //mObject.startDir.print();
                m4d::vec4 oldDir = m4d::vec4(0.0,mObject.startDir[0],mObject.startDir[1],mObject.startDir[2]);
                //oldDir.printS();
                m4d::vec4 newDir;
                m4d::TransCoordinates::transCoordCartSph(oldDir,newDir);
                //newDir.printS();
                mObject.ksi = newDir[3] * RAD_TO_DEG;
                mObject.chi = newDir[2] * RAD_TO_DEG;
                //std::cerr << mObject.ksi << " " << mObject.chi << std::endl;
            }

            if (getfield(L,"type",pVal)) {
                mObject.type = (m4d::enum_geodesic_type)((int)pVal);
                mObject.geodSolver->setGeodesicType(mObject.type);
            }

            if (getfield(L,"type",sVal)) {
                if (sVal.compare("lightlike")==0 || sVal.compare("null")==0) {
                    mObject.type = m4d::enum_geodesic_lightlike;
                }
                else if (sVal.compare("spacelike")==0) {
                    mObject.type = m4d::enum_geodesic_spacelike;
                }
                else if (sVal.compare("timelike")==0) {
                    mObject.type = m4d::enum_geodesic_timelike;
                }
                else if (sVal.compare("lightlikeSachs")==0) {
                    mObject.type = m4d::enum_geodesic_lightlike_sachs;
                }
                mObject.geodSolver->setGeodesicType(mObject.type);
            }

            if (getfield(L,"timedir",pVal)) {
                mObject.timeDirection = (int)pVal;
            }

            if (getfield(L,"timedir",sVal)) {
                if (sVal.compare("forward")==0) {
                    mObject.timeDirection = 1.0;
                }
                else if (sVal.compare("backward")==0) {
                    mObject.timeDirection = -1.0;
                }
                else {
                    mObject.timeDirection = 0.0;
                }
            }

            if (getfield(L,"vel",pVal)) {
                mObject.vel = pVal;
            }

            if (getfield(L,"maxpoints",pVal)) {
                mObject.maxNumPoints = (unsigned int)pVal;
            }

            //mObject.geodSolver->calculateGeodesic(mObject.startPos
            double eta = 1.0;
            double y0  = 1.0;

            if (mObject.type==m4d::enum_geodesic_timelike && fabs(mObject.vel)<1.0) {
                double beta = mObject.vel/mObject.currMetric->speed_of_light();
                double gm = 1.0-beta*beta;
                if (gm>0.0) {
                    y0  = mObject.currMetric->speed_of_light()/sqrt(gm);
                    eta = beta*y0;
                }
            }

            m4d::vec4 b[4];
            for(int i=0; i<4; i++) {
                for(int k=0; k<4; k++) {
                    b[i] += mObject.lorentz.getElem(i,k)*mObject.base[k];
                }
            }

            if (mObject.type == m4d::enum_geodesic_spacelike) {
                mObject.timeDirection = 0.0;
            }

            // Initial direction with respect to natural local tetrad:
            //m4d::vec4 locDir  = DEF_LUA_SIGNUM(mObject.timeDirection)*y0*b[0] + eta*(mObject.startDir[0]*b[1] + mObject.startDir[1]*b[2] + mObject.startDir[2]*b[3]);
            m4d::vec4 locDir  = mObject.timeDirection*y0*b[0] + eta*(mObject.startDir[0]*b[1] + mObject.startDir[1]*b[2] + mObject.startDir[2]*b[3]);
            m4d::vec3 nullDir = m4d::vec3(mObject.startDir[0], mObject.startDir[1], mObject.startDir[2]);
            m4d::vec4 coDir;
            mObject.currMetric->localToCoord( mObject.startPos, locDir, coDir, mObject.tetradType );
            mObject.geodSolver->setAffineParamStep( mObject.stepsize );
            mObject.coordDir = coDir;

            m4d::vec3 locX = m4d::vec3(1.0,0.0,0.0);
            m4d::vec3 locY = m4d::vec3(0.0,1.0,0.0);
            m4d::vec3 locZ = m4d::vec3(0.0,0.0,1.0);
            m4d::vec3 lX,lY,lZ;
            lX = locX; lY = locY; lZ = locZ;

            //m4d::enum_break_condition  breakCond = m4d::enum_break_none;


            if (mObject.type == m4d::enum_geodesic_lightlike ||
                    mObject.type == m4d::enum_geodesic_timelike ||
                    mObject.type == m4d::enum_geodesic_spacelike) {
                mObject.geodSolver->calculateGeodesic(mObject.startPos,coDir,mObject.maxNumPoints,mObject.points,mObject.dirs,mObject.lambda);
            }
            else if (mObject.type == m4d::enum_geodesic_lightlike_sachs) {
                mObject.geodSolver->calcSachsJacobi(mObject.startPos,coDir,nullDir,lX,lY,lZ,b[0],b[1],b[2],b[3],mObject.tetradType,mObject.maxNumPoints,mObject.points,mObject.dirs,mObject.lambda,mObject.sachs1,mObject.sachs2,mObject.jacobi,mObject.maxJacobi);
            }

            //fprintf(stderr,"%s\n",m4d::stl_break_condition[breakCond]);
        }
    }
    lua_pushnumber(L,mObject.points.size());
    return 1;
}


/**
 * @brief Interpolate position by coordinates.
 *
 *   After calculating a geodesic, a position along it can be interpolated by means of
 *   a coordinate value. For that, two values, the coordinate index (0,1,2,3) and the
 *   coordinate value has to be given as parameters.
 * @param L Pointer to lua state.
 * @return 0 : if error occurred.\n
 *         1 : the interpolated 4-point is put onto the lua stack.
 */
int interpolatePosByCoord( lua_State *L ) {
    int numParams = lua_gettop(L);
    if (numParams!=2 || mObject.lambda.size()==0 || mObject.points.size()!=mObject.lambda.size()) {
        // TODO: error handling
        return 0;
    }

    int coordNum = 0;
    double s = 0.0;

    if (!lua_isnil(L,-2)) {
        if (lua_isinteger(L,-2)) {
            coordNum = lua_tointeger(L,-2);

            if (!lua_isnil(L,-1)) {
                if (lua_isnumber(L,-1)) {
                    s = lua_tonumber(L,-1);

                    unsigned int j=0;
                    double l1,l2,x;
                    m4d::vec4 pos = m4d::vec4();

                    while(j < mObject.points.size()-1) {
                        l1 = mObject.points[j][coordNum];
                        l2 = mObject.points[j+1][coordNum];

                        if ((l1-s)*(l2-s) < 0.0) {
                            x = (s - l1)/(l2 - l1);
                            pos = mObject.points[j]*(1.0-x) + mObject.points[j+1]*x;

                            lua_newtable(L);
                            for(int i=1; i<=4; i++) {
                                lua_pushnumber(L,i);
                                lua_pushnumber(L,pos[i-1]);
                                lua_settable(L,-3);
                            }
                            return 1;
                        }
                        j++;
                    }
                }
            }
        }
    }
    return 0;
}


/**
 * @brief Interpolate position by affine parameter.
 *
 *   After calculating a geodesic, a position along it can be interpolated by means of
 *   the affine parameter which has to be given as parameter.
 * @param L Pointer to lua state.
 * @return 0 : if error occurred.\n
 *         1 : the interpolated 4-point is put onto the lua stack.
 */
int interpolatePosByLambda( lua_State *L ) {
    int numParams = lua_gettop(L);
    if (numParams==0 || mObject.lambda.size()==0 || mObject.points.size()!=mObject.lambda.size()) {
        // TODO: error handling
        return 0;
    }

    double lambda = 0.0;
    if (lua_isnil(L,-1) || !lua_isnumber(L,-1)) {
        mlua_error(L,"InterpolatePosByLambda: argument must be a number.\n");
        return 0;
    }
    lambda = lua_tonumber(L,-1);

    unsigned int j=0;
    double l1,l2,x;
    m4d::vec4 pos = m4d::vec4();

    while(j < mObject.lambda.size()-1) {
        l1 = mObject.lambda[j];
        l2 = mObject.lambda[j+1];
        if ((l1-lambda)*(l2-lambda) < 0.0) {
            x = (lambda - l1)/(l2 - l1);
            pos = mObject.points[j]*(1.0-x) + mObject.points[j+1]*x;

            lua_newtable(L);
            for(int i=1; i<=4; i++) {
                lua_pushnumber(L,i);
                lua_pushnumber(L,pos[i-1]);
                lua_settable(L,-3);
            }
            return 1;
        }
        j++;
    }

    return 0;
}


/**
 * @brief Get current affine parameter value.
 * @param L Pointer to lua state.
 * @return 0 : if error occurred.\n
 *         1 : the affine parameter is put onto the lua stack.
 */
int getCurrentLambda(lua_State *L) {
    if (mObject.geodSolver == NULL) {
        fprintf(stderr,"No Geodesic solver initialized!\n");
        return 0;
    }

    double lambda = mObject.geodSolver->getAffineParam();
    lua_pushnumber(L,lambda);
    return 1;
}


/**
 * @brief getCurrentPos
 * @param L
 * @return
 */
int getCurrentPos(lua_State* L) {
    if (mObject.geodSolver == NULL) {
        fprintf(stderr,"No Geodesic solver initialized!\n");
        return 0;
    }

    double pos[4];
    mObject.geodSolver->getPosition(pos);

    /*
    lua_newtable(L);
    for(int i=1; i<=4; i++) {
        lua_pushnumber(L,i);
        lua_pushnumber(L,pos[i-1]);
        lua_settable(L,-3);
    }
    return 1;
    */
    return pushVec(L,pos);
}


/**
 * @brief getCurrentDir
 * @param L
 * @return
 */
int getCurrentDir(lua_State* L) {
    if (mObject.geodSolver == NULL) {
        fprintf(stderr,"No Geodesic solver initialized!\n");
        return 0;
    }

    double dir[4];
    mObject.geodSolver->getDirection(dir);

    /*
    lua_newtable(L);
    for(int i=1; i<=4; i++) {
        lua_pushnumber(L,i);
        lua_pushnumber(L,dir[i-1]);
        lua_settable(L,-3);
    }
    return 1;
    */
    return pushVec(L,dir);
}


/**
 * @brief Save geodesic to file.
 * @param L
 * @return
 */
int saveGeodesic(lua_State* L) {
    int numParams = lua_gettop(L);
    if (numParams==0) {
        return 0;
    }

    std::string filename;
    if (!lua_isnil(L,-1)) {
        if (lua_isstring(L,-1)) {
            filename = lua_tostring(L,-1);
            fprintf(stderr,"file: %s\n",filename.c_str());
            FILE* fptr;
            if ((fptr=fopen(filename.c_str(),"w"))==NULL) {
                return 0;
            }

            for(unsigned int i=0; i<mObject.points.size(); i++) {
                fprintf(fptr,"%5d %12.8f ",i,mObject.lambda[i]);mObject.points[i].printS(fptr);
            }
            fclose(fptr);
        }
    }
    return 0;
}
