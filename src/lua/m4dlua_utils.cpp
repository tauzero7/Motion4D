/* -------------------------------------------------------------------------------
    m4dlua_utils.cpp

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

#include "lua/m4dlua_utils.h"
#include "m4dGlobalDefs.h"

void lua_reg_utils(lua_State *L) {
    lua_register(L, "PrintObjectSettings", printObjectSettings);

    lua_register(L, "GetCameraView", getCameraView);

    lua_register(L, "GetGeodPoint", getGeodPoint);
    lua_register(L, "GetGeodDir", getGeodDir);
    lua_register(L, "GetGeodLambda", getGeodLambda);
    lua_register(L, "GetNumPoints", getNumPoints);

    lua_register(L, "LinInterpolVec", linInterpolVec);

    lua_register(L, "PrintVec", printVec);
    lua_register(L, "PrintMat", printMat);
}

/**
 * @brief stackDump
 * @param L
 */
void stackDump(lua_State *L)  {
    int i;
    int top = lua_gettop(L);
    for (i = 1; i <= top; i++) {  /* repeat for each level */
        int t = lua_type(L, i);
        switch (t) {

            case LUA_TSTRING: { /* strings */
                printf("`%s'", lua_tostring(L, i));
                break;
            }
            case LUA_TBOOLEAN: { /* booleans */
                printf(lua_toboolean(L, i) ? "true" : "false");
                break;
            }
            case LUA_TNUMBER: { /* numbers */
                printf("%g", lua_tonumber(L, i));
                break;
            }
            case LUA_TTABLE: {
                printf("table");
                break;
            }
            default: { /* other values */
                printf("%s", lua_typename(L, t));
                break;
            }

        }
        printf("  ");  /* put a separator */
    }
    printf("\n");  /* end the listing */
}

/**
 * @brief error
 * @param L
 * @param fmt
 */
void mlua_error(lua_State *L, const char *fmt, ...) {
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr,fmt,argp);
    va_end(argp);
    lua_close(L);
    //exit(EXIT_FAILURE);
}


/**
 * @brief getTableValues
 * @param L
 * @param vec
 */
bool getTableValues(lua_State *L, std::vector<double> &vec, int rel) {
    vec.clear();

    int numParams = lua_gettop(L);
    if (numParams==0 || -rel > numParams) {
        return false;
    }

    if (!lua_isnil(L,rel)) {
        if (lua_istable(L,rel)) {
            for(unsigned int i=1; i<=lua_rawlen(L,rel); i++) {
                lua_rawgeti(L,rel,i);
                int top = lua_gettop(L);
                double v = lua_tonumber(L,top);
                vec.push_back(v);
                lua_pop(L,1);
            }
            return true;
        }
    }
    return false;
}



bool getfield( lua_State *L, const char *key, int &val ) {
    bool isValid = true;
    lua_pushstring(L, key);
    lua_gettable(L, -2);

    //std::cerr << key << std::endl;
    if (!lua_isinteger(L, -1)) {
        isValid = false;
    }
    else {
        val = (int)lua_tointeger(L, -1);
    }
    lua_pop(L, 1);  /* remove number */
    return isValid;
}


/**
 * @brief getfield
 * @param L
 * @param key
 * @param val
 * @return
 */
bool  getfield ( lua_State *L, const char *key, double &val ) {
    bool isValid = true;
    lua_pushstring(L, key);
    lua_gettable(L, -2);

    //std::cerr << key << std::endl;
    if (!lua_isnumber(L, -1)) {
        isValid = false;
    }
    else {
        val = (double)lua_tonumber(L, -1);
    }
    lua_pop(L, 1);  /* remove number */
    return isValid;
}


/**
 * @brief getfield
 * @param L
 * @param key
 * @param str
 * @return
 */
bool getfield(lua_State *L, const char *key, std::string &str) {
    bool isValid = true;
    lua_pushstring(L, key);
    lua_gettable(L, -2);
    if (!lua_isstring(L,-1)) {
        isValid = false;
    }
    else {
        str = lua_tostring(L,-1);
    }
    lua_pop(L, 1);
    return isValid;
}


/**
 * @brief getfield
 * @param L
 * @param key
 * @param vals
 * @return
 */
bool getfield(lua_State *L, const char *key, std::vector<double> &vals) {
    bool isValid = true;
    lua_pushstring(L, key);
    lua_gettable(L, -2);
    if (!lua_istable(L,-1)) {
        isValid = false;
    }
    else {
        for(unsigned int i=1; i<=lua_rawlen(L,-1); i++) {
            lua_rawgeti(L,-1,i);
            int top = lua_gettop(L);
            double v = lua_tonumber(L,top);
            vals.push_back(v);
            lua_pop(L,1);
        }
    }
    lua_pop(L, 1);
    return isValid;
}


/**
 * @brief getvector
 * @param L
 * @param name
 * @param vals
 * @return
 */
int getvector(lua_State* L, std::string name, std::vector<double> &vals) {
    if (!vals.empty()) {
        vals.clear();
    }

    lua_getglobal( L, name.c_str() );
    if (!lua_istable(L,-1)) {
        return 0;
    }

    for(unsigned int i=1; i<=lua_rawlen(L,-1); i++) {
        lua_rawgeti(L,-1,i);
        int top = lua_gettop(L);
        double v = lua_tonumber(L,top);
        vals.push_back(v);
        lua_pop(L,1);
    }
    return vals.size();
}


/**
 * @brief Push vector to lua stack.
 *
 *   The double array must be of size 4!
 * @param L
 * @param vec
 * @return
 */
int pushVec(lua_State *L, double* vec) {
    lua_newtable(L);
    for(int i=1; i<=4; i++) {
        lua_pushnumber(L,i);
        lua_pushnumber(L,vec[i-1]);
        lua_settable(L,-3);
    }
    return 1;
}


int pushVec(lua_State *L, m4d::vec3 &vec) {
    lua_newtable(L);
    for(int i=1; i<=3; i++) {
        lua_pushnumber(L,i);
        lua_pushnumber(L,vec[i-1]);
        lua_settable(L,-2);
    }
    return 1;
}


int pushVec(lua_State *L, m4d::vec4 &vec) {    
    lua_newtable(L);
    for(int i=1; i<=4; i++) {
        lua_pushnumber(L,i);
        lua_pushnumber(L,vec[i-1]);
        lua_settable(L,-3);
    }
    return 1;
}


int pushMat(lua_State *L, m4d::mat4 &mat) {
    lua_createtable(L,4,0);
    for(int y=1; y<=4; y++) {
        lua_pushnumber(L,y);
        lua_createtable(L,0,4);
        for(int x=1; x<=4; x++) {
            lua_pushnumber(L,x);
            lua_pushnumber(L,mat.getElem(y-1,x-1));
            lua_settable(L,-3);
        }
        lua_settable(L,-3);
    }
    return 1;
}



int getCameraView(lua_State *L) {
    if (lua_isnil(L,-1) || !lua_istable(L,-1)) {
        mlua_error(L,"GetCameraView....\n");
        return 0;
    }


    int resX, resY;
    double fovX, fovY;
    bool haveResX = getfield(L, "resX", resX);
    bool haveResY = getfield(L, "resY", resY);
    bool haveFovX = getfield(L, "fovX", fovX);
    bool haveFovY = getfield(L, "fovY", fovY);

    if (haveFovX && haveFovY) {
        if (haveResX && haveResY) {
            fprintf(stderr,"GetCameraView: one entry should be left.\n");
            return 0;
        }
        else if (!haveResX && !haveResY) {
            fprintf(stderr,"GetCameraView: either resX or resY have to be defined.\n");
            return 0;
        }
    }


    if (haveResX && haveResY) {
        if (haveFovX) {
            fovY = m4d::degree(2.0 * atan(resY/static_cast<double>(resX) * tan(m4d::radians(fovX*0.5))));
            lua_pushnumber(L,fovY);
            return 1;
        }
        else {
            fovX = m4d::degree(2.0 * atan(resX/static_cast<double>(resY) * tan(m4d::radians(fovY*0.5))));
            lua_pushnumber(L,fovX);
            return 1;
        }
    }

    if (haveResX && haveFovX && haveFovY) {
        resY = static_cast<int>(resX * tan(m4d::radians(fovY*0.5)) / tan(m4d::radians(fovX*0.5)));
        lua_pushinteger(L,resY);
        return 1;
    }

    if (haveResY && haveFovX && haveFovY) {
        resX = static_cast<int>(resY * tan(m4d::radians(fovX*0.5)) / tan(m4d::radians(fovY*0.5)));
        lua_pushinteger(L,resX);
        return 1;
    }

    fprintf(stderr,"GetCameraView: missing arguments.\n");
    return 0;
}


/**
 * @brief getGeodPoint
 * @param L
 * @return
 */
int getGeodPoint(lua_State *L) {
    if (lua_isnil(L,-1) || !lua_isinteger(L,-1)) {
        mlua_error(L,"GetGeodPoint needs valid point index.\n");
        return 0;
    }

    int idx = lua_tointeger(L,-1);
    if (idx<0 || idx >= (int)mObject.points.size()) {
        mlua_error(L,"GetGeodPoint needs valid point index.\n");
        return 0;
    }

/*
    lua_newtable(L);
    for(int i=1; i<=4; i++) {
        lua_pushnumber(L,i);
        lua_pushnumber(L,mObject.points[idx][i-1]);
        lua_settable(L,-3);
    }
    return 1;
*/
    return pushVec(L,mObject.points[idx]);
}


/**
 * @brief getGeodDir
 * @param L
 * @return
 */
int getGeodDir(lua_State *L) {
    if (lua_isnil(L,-1) || !lua_isinteger(L,-1)) {
        mlua_error(L,"GetGeodPoint needs valid point index.\n");
        return 0;
    }

    int idx = lua_tointeger(L,-1);
    if (idx<0 || idx >= (int)mObject.dirs.size()) {
        mlua_error(L,"GetGeodPoint needs valid point index.\n");
        return 0;
    }

    /*
    lua_newtable(L);
    for(int i=1; i<=4; i++) {
        lua_pushnumber(L,i);
        lua_pushnumber(L,mObject.dirs[idx][i-1]);
        lua_settable(L,-3);
    }
    return 1;
    */
    return pushVec(L,mObject.dirs[idx]);
}


/**
 * @brief getGeodLambda
 * @param L
 * @return
 */
int getGeodLambda(lua_State *L) {
    if (lua_isnil(L,-1) || !lua_isinteger(L,-1)) {
        mlua_error(L,"GetGeodPoint needs valid point index.\n");
        return 0;
    }

    int idx = lua_tointeger(L,-1);
    if (idx<0 || idx >= (int)mObject.points.size()) {
        mlua_error(L,"GetGeodPoint needs valid point index.\n");
        return 0;
    }

    lua_pushnumber(L,mObject.lambda[idx]);
    return 1;
}


/**
 * @brief getNumPoints
 * @param L
 * @return
 */
int getNumPoints(lua_State *L) {
    lua_pushinteger(L,mObject.points.size());
    return 1;
}


/**
 * @brief linInterpolVec
 * @param L
 * @return
 */
int linInterpolVec(lua_State *L) {
    int numParams = lua_gettop(L);
    if (numParams != 3) {
        fprintf(stderr,"LinInterpolVec: missing arguments\n");
        return 0;
    }

    std::vector<double> v1, v2;
    double x;

    if (lua_isnil(L,-3) || !lua_istable(L,-3) ||
            lua_isnil(L,-2) || !lua_istable(L,-2) ||
            lua_isnil(L,-1) || !lua_isnumber(L,-1)) {
        fprintf(stderr,"LinInterpolVec: wrong arguments\n");
        return 0;
    }

    getTableValues(L,v1,-3);
    getTableValues(L,v2,-2);
    x = lua_tonumber(L,-1);

    if ((v1.size() != v2.size()) || v1.size() != 4) {
        fprintf(stderr,"LinInterpolVec: Vectors have to be of size 4!\n");
        return 0;
    }

    m4d::vec4 p1 = m4d::vec4(v1[0],v1[1],v1[2],v1[3]);
    m4d::vec4 p2 = m4d::vec4(v2[0],v2[1],v2[2],v2[3]);
    m4d::vec4 ip = (1.0-x) * p1 + x * p2;

    return pushVec(L,ip);
}


/**
 * @brief printObjectSettings
 * @return
 */
int printObjectSettings(lua_State *) {
    mObject.printSettings(stdout);
    return 0;
}


/**
 * @brief Print vector
 *      Use PrintVec(vec, "%12.8f ").
 *      The second argument is optional.
 *
 * @param L
 * @return
 */
int printVec(lua_State *L) {
    std::string format = "%8.4f ";
    int numParams = lua_gettop(L);
    int relpos = -1;
    if (numParams==2) {
        if (lua_isstring(L,-1)) {
            format = lua_tostring(L,-1);
            relpos = -2;
        }
    }

    std::vector<double> vec;
    getTableValues(L,vec,relpos);
    if (vec.empty()) {
        return 0;
    }

    for(unsigned int i=0; i < vec.size(); i++) {
        printf(format.c_str(),vec[i]);
    }
    printf("\n");
    return 0;
}


/**
 * @brief printMat
 * @param L
 * @return
 */
int printMat(lua_State *L) {
    std::string format = "%8.4f ";
    int numParams = lua_gettop(L);
    int relpos = -1;
    if (numParams==2) {
        if (lua_isstring(L,-1)) {
            format = lua_tostring(L,-1);
            relpos = -2;
        }
    }

    std::vector<std::vector<double> > rows;

    if (!lua_isnil(L,relpos) && lua_istable(L,relpos)) {
        for(unsigned int i=1; i<=lua_rawlen(L,relpos); i++) {
            std::vector<double> vals;
            lua_pushnumber(L, i);
            lua_gettable(L, relpos-1);
            if (lua_istable(L,-1)) {
                for(unsigned int j=1; j<=lua_rawlen(L,-1); j++) {
                    lua_rawgeti(L,-1,j);
                    int top = lua_gettop(L);
                    double v = lua_tonumber(L,top);
                    vals.push_back(v);
                    lua_pop(L,1);
                }
            }
            lua_pop(L, 1);
            rows.push_back(vals);
        }
    }

    for(unsigned int i=0; i<rows.size(); i++) {
        for(unsigned int j=0; j<rows[i].size(); j++) {
            printf(format.c_str(),rows[i][j]);
        }
        printf("\n");
    }

    return 0;
}
