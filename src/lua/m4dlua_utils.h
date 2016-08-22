/* -------------------------------------------------------------------------------
    m4dlua_utils.h

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

#ifndef M4D_LUA_UTILS_H
#define M4D_LUA_UTILS_H

#include "lua/m4dlua.h"

// ------------------------------------------------------------------------------
//  internal functions
// ------------------------------------------------------------------------------

void  M4D_CALL lua_reg_utils(lua_State *L);

void  stackDump( lua_State *L);
void  mlua_error( lua_State *L, const char *fmt, ...);

bool  getTableValues( lua_State *L, std::vector<double> &vec, int rel = -1 );

bool  getfield  ( lua_State *L, const char *key, int &val );
bool  getfield  ( lua_State *L, const char *key, double &val );
bool  getfield  ( lua_State *L, const char *key, std::string &str );
bool  getfield  ( lua_State *L, const char *key, std::vector<double> &vals );
int   getvector ( lua_State* L, std::string name, std::vector<double> &vals );

int   pushVec(lua_State *L, double* vec);
int   pushVec(lua_State *L, m4d::vec3 &vec);
int   pushVec(lua_State *L, m4d::vec4 &vec);

int   pushMat(lua_State *L, m4d::mat4 &mat);

// ------------------------------------------------------------------------------
//  functions available in lua
// ------------------------------------------------------------------------------
int   getCameraView(lua_State *L);

int   getGeodPoint(lua_State *L);
int   getGeodDir(lua_State *L);
int   getGeodLambda(lua_State *L);
int   getTrajE0(lua_State *L);
int   getTrajE1(lua_State *L);
int   getTrajE2(lua_State *L);
int   getTrajE3(lua_State *L);
int   getTrajEi(lua_State *L, unsigned int i);

int   getNumPoints(lua_State *L);

int   linInterpolVec(lua_State *L);

int   printObjectSettings(lua_State *L);
int   printVec(lua_State *L);
int   printMat(lua_State *L);


#endif // M4D_LUA_UTILS_H
