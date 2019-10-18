/* -------------------------------------------------------------------------------
    m4dlua_metric.cpp

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

#include "lua/m4dlua_metric.h"

#include "m4dGlobalDefs.h"
#include "math/TransCoordinates.h"

extern m4d::Object  mObject;

/**
 * @brief lua_reg_metric
 * @param L
 */
void lua_reg_metric(lua_State* L) {
    lua_register(L, "PrintMetricDB",   printMetricDB);

    lua_register(L, "SetMetric",       setMetric);
    lua_register(L, "SetMetricParam",  setMetricParam);
    lua_register(L, "SetMetricParams", setMetricParams);
    lua_register(L, "GetMetricParam",  getMetricParam);
    lua_register(L, "PrintMetric",     printMetric);

    lua_register(L, "LocalToCoord", localToCoord);
    lua_register(L, "CoordToLocal", coordToLocal);

    lua_register(L, "ToPseudoCoords", transToPseudoCoords);
    lua_register(L, "CoordTrans", coordTrans);

    lua_register(L, "CalcProduct", calcProduct);
    lua_register(L, "CalcTidalMatrix", calcTidalMatrix);
}


int printMetricDB(lua_State* ) {
    mObject.metricDB.printMetricList();
    return 0;
}


int setMetric(lua_State* L) {
    int numParams = lua_gettop(L);
    if (numParams==0) {
        mlua_error(L,"SetMetric expects at least a metric name.\n");
        return 0;
    }

    if (!lua_isnil(L,-numParams)) {
        if (lua_isstring(L,-numParams)) {
            std::string metricName = std::string(lua_tostring(L,-numParams));
            if (!mObject.setMetric(metricName.c_str())) {
                mlua_error(L,"SetMetric: no valid metric set.\n");
                return 0;
            }
        }
    } else {
        mlua_error(L,"SetMetric expects a metric name as first argument.\n");
        return 0;
    }

    if (numParams==2 && !lua_isnil(L,-1)) {
        if (lua_istable(L,-numParams+1)) {

            std::vector<std::string> paramNames;
            mObject.currMetric->getParamNames(paramNames);
            double paramVal;

            for(unsigned int i=0; i<paramNames.size(); i++) {
                if (getfield(L,paramNames[i].c_str(),paramVal))
                    mObject.currMetric->setParam(paramNames[i].c_str(),paramVal);
            }
        }
    }

    if (numParams>2){
        mlua_error(L,"SetMetric: Wrong number of parameters!\n");
        lua_pop(L,lua_gettop(L));
    }
    return 0;
}


int setMetricParam(lua_State* L) {
    int numParams = lua_gettop(L);
    if (numParams!=2) {
        lua_pop(L,numParams);
        return 0;
    }

    std::string paramName;
    double      paramValue = 0.0;

    if (!lua_isnil(L,-2)) {
        if (lua_isstring(L,-2)) {
            paramName = std::string(lua_tostring(L,-2));
        }
    }

    if (!lua_isnil(L,-1)) {
        if (lua_isnumber(L,-1)) {
            paramValue = lua_tonumber(L,-1);
        }
    }

    if (mObject.currMetric!=NULL) {
        mObject.currMetric->setParam(paramName.c_str(),paramValue);
    }
    return 0;
}


/**
 * @brief setMetricParams
 * @param L
 * @return
 */
int setMetricParams(lua_State* L) {
    if (mObject.currMetric == NULL) {
        // TODO: error handling
        return 0;
    }

    std::vector<std::string> paramNames;
    mObject.currMetric->getParamNames(paramNames);
    double paramVal;

    if (!lua_isnil(L,-1)) {
        if (lua_istable(L,-1)) {
            for(unsigned int i=0; i<paramNames.size(); i++) {
                if (getfield(L,paramNames[i].c_str(),paramVal))
                    mObject.currMetric->setParam(paramNames[i].c_str(),paramVal);
            }
        }
    }
    return 0;
}


int printMetric(lua_State* ) {
    if (mObject.currMetric!=NULL) {
        mObject.currMetric->printF(stdout);
    }
    return 0;
}


int getMetricParam(lua_State *L) {
    if (mObject.currMetric == NULL) {
        // TODO: error handling
        return 0;
    }

    if (!lua_isnil(L,-1)) {
        if (lua_isstring(L,-1)) {
            std::string paramName = lua_tostring(L,-1);
            double paramValue;
            if (mObject.currMetric->getParam(paramName.c_str(),paramValue)) {
                lua_pushnumber(L,paramValue);
                return 1;
            }
        }
    }
    return 0;
}


int localToCoord(lua_State *L) {
    if (mObject.currMetric == NULL) {
        fprintf(stderr,"LocalToCoord: No metric set!\n");
        return 0;
    }

    if (lua_isnil(L,-1) || !lua_istable(L,-1)) {
        return 0;
    }

    std::vector<double> dpos;
    if (!getfield(L, "pos", dpos) || dpos.size()<4) {
        fprintf(stderr,"LocalToCoord: pos parameter is missing!\n");
        return 0;
    }
    m4d::vec4 pos = m4d::vec4(dpos[0],dpos[1],dpos[2],dpos[3]);

    std::vector<double> dvec;
    if (!getfield(L, "vec", dvec) || dvec.size()<4) {
        fprintf(stderr,"LocalToCoord: vec parameter is missing!\n");
        return 0;
    }
    m4d::vec4 vec = m4d::vec4(dvec[0],dvec[1],dvec[2],dvec[3]);

    m4d::vec4 nvec;
    std::string sType;
    m4d::enum_nat_tetrad_type eType = m4d::enum_nat_tetrad_default;
    if (getfield(L, "type", sType)) {
        int nType = m4d::find_nat_tetrad_type(sType.c_str());
        if (nType>=0) {
            eType = static_cast<m4d::enum_nat_tetrad_type>(nType);
        }
    }

    mObject.currMetric->localToCoord(pos,vec,nvec,eType);
    return pushVec(L,nvec);
}


int  coordToLocal(lua_State *L) {
    if (mObject.currMetric == NULL) {
        fprintf(stderr,"CoordToLocal: No metric set!");
        return 0;
    }

    if (lua_isnil(L,-1) || !lua_istable(L,-1)) {
        return 0;
    }

    std::vector<double> dpos;
    if (!getfield(L, "pos", dpos) || dpos.size()<4) {
        // error
        return 0;
    }
    m4d::vec4 pos = m4d::vec4(dpos[0],dpos[1],dpos[2],dpos[3]);

    std::vector<double> dvec;
    if (!getfield(L, "vec", dvec) || dvec.size()<4) {
        fprintf(stderr,"LocalToCoord: vec parameter is missing!\n");
        return 0;
    }
    m4d::vec4 vec = m4d::vec4(dvec[0],dvec[1],dvec[2],dvec[3]);

    m4d::vec4 nvec;
    std::string sType;
    m4d::enum_nat_tetrad_type eType = m4d::enum_nat_tetrad_default;
    if (getfield(L, "type", sType)) {
        int nType = m4d::find_nat_tetrad_type(sType.c_str());
        if (nType>=0) {
            eType = static_cast<m4d::enum_nat_tetrad_type>(nType);
        }
    }

    mObject.currMetric->coordToLocal(pos,vec,nvec,eType);
    return pushVec(L,nvec);
}


/**
 * @brief transToPseudoCoords
 * @param L
 * @return
 */
int transToPseudoCoords( lua_State *L ) {
    if (mObject.currMetric == NULL) {
        // TODO: error handling
        return 0;
    }

    std::vector<double> vec;
    if (getTableValues(L,vec)) {
        if (vec.size()==4) {
            m4d::vec4 pos = m4d::vec4(vec[0], vec[1], vec[2] , vec[3]);
            m4d::vec4 cpos;
            mObject.currMetric->transToPseudoCart(pos,cpos);

            lua_newtable(L);
            for(int i=1; i<=4; i++) {
                lua_pushnumber(L,i);
                lua_pushnumber(L,cpos[i-1]);
                lua_settable(L,-3);
            }
            return 1;
        }
    }
    return 0;
}


int coordTrans(lua_State *L) {
    int numParams = lua_gettop(L);
    if (numParams != 3 || !lua_istable(L,-3) || !lua_isstring(L,-2) || !lua_isstring(L,-1)) {
        mlua_error(L, "CoordTrans needs 4-vector, fromCoord, toCoord as arguments.\n");
        return 0;
    }

    std::string fromCoordName = lua_tostring(L,-2);
    std::string toCoordName = lua_tostring(L,-1);

    m4d::enum_coordinate_type fromCoordType = m4d::enum_coordinate_cartesian;
    m4d::enum_coordinate_type toCoordType = m4d::enum_coordinate_cartesian;

    unsigned int i=0;
    while(i < m4d::NUM_ENUM_COORDINATE_TYPES) {
        if (fromCoordName.compare(m4d::stl_coordinate_types[i])==0) {
            fromCoordType = (m4d::enum_coordinate_type)i;
            break;
        }
        i++;
    }

    i=0;
    while(i < m4d::NUM_ENUM_COORDINATE_TYPES) {
        if (toCoordName.compare(m4d::stl_coordinate_types[i])==0) {
            toCoordType = (m4d::enum_coordinate_type)i;
            break;
        }
        i++;
    }

    std::vector<double> vec;
    if (getTableValues(L,vec,-3)) {
        if (vec.size()==4) {
            m4d::vec4 oldVec = m4d::vec4(vec[0],vec[1],vec[2],vec[3]);
            m4d::vec4 newVec;
            m4d::TransCoordinates::coordTransf(fromCoordType,oldVec,toCoordType,newVec);

            lua_newtable(L);
            for(int i=1; i<=4; i++) {
                lua_pushnumber(L,i);
                lua_pushnumber(L,newVec[i-1]);
                lua_settable(L,-3);
            }
            return 1;
        }
    }
    return 0;
}


int calcProduct( lua_State *L ) {
    if (mObject.currMetric == NULL) {
        fprintf(stderr,"CalcProduct: metric not set!\n");
        return 0;
    }

    if (lua_isnil(L,-1) || !lua_istable(L,-1)) {
        fprintf(stderr, "CalcProduct: needs parameters\n");
        return 0;
    }

    std::vector<double>  dpos;
    if (!getfield(L, "pos", dpos) || dpos.size()<4) {
        fprintf(stderr, "CalcProduct: position vector is missing\n");
        return 0;
    }

    std::vector<double>  du;
    if (!getfield(L, "u", du) || du.size()<4) {
        fprintf(stderr, "CalcProduct: vector u is missing\n");
        return 0;
    }

    std::vector<double>  dv;
    if (!getfield(L, "v", dv) || dv.size()<4) {
        fprintf(stderr, "CalcProduct: vector v is missing\n");
        return 0;
    }

    m4d::vec4 pos = m4d::vec4(&dpos[0], 4);
    m4d::vec4 u = m4d::vec4(&du[0], 4);
    m4d::vec4 v = m4d::vec4(&dv[0], 4);

    double prod;
    mObject.currMetric->calcProduct(pos, u, v, prod);

    lua_pushnumber(L, prod);
    return 1;
}


int calcTidalMatrix( lua_State *L ) {
    if (mObject.currMetric == NULL) {
        fprintf(stderr,"CalcTidalMatrix: metric not set!\n");
        return 0;
    }


    if (lua_isnil(L,-1) || !lua_istable(L,-1)) {
        fprintf(stderr,"CalcTidalMatrix: at least 'pos' parameter has to be set!\n");
        return 0;
    }

    std::vector<double> dpos;
    if (!getfield(L, "pos", dpos) || dpos.size()<4) {
        fprintf(stderr,"CalcTidalMatrix: pos parameter is missing!\n");
        return 0;
    }
    m4d::vec4 pos = m4d::vec4(dpos[0],dpos[1],dpos[2],dpos[3]);

    std::string sType;
    m4d::enum_nat_tetrad_type eType = m4d::enum_nat_tetrad_default;
    if (getfield(L, "type", sType)) {
        int nType = m4d::find_nat_tetrad_type(sType.c_str());
        if (nType>=0) {
            eType = static_cast<m4d::enum_nat_tetrad_type>(nType);
        }
    }

    m4d::vec4 e0,e1,e2,e3;
    mObject.currMetric->getNatTetrad(pos,e0,e1,e2,e3,eType);
    mObject.currMetric->calcTidalMatrix(pos,e0,e1,e2,e3);

    //double mm[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    m4d::mat4 tm;
    mObject.currMetric->getTidalMatrix(tm);

    return pushMat(L,tm);
}
