/* -------------------------------------------------------------------------------
    m4dlua_solver.h

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

#ifndef M4D_LUA_SOLVER_H
#define M4D_LUA_SOLVER_H

#include "lua/m4dlua.h"
#include "lua/m4dlua_utils.h"

void  lua_reg_solver(lua_State *L);

int  printSolverDB  ( lua_State *L );
int  printGeodTypes ( lua_State *L );

int  setGeodSolver ( lua_State *L );
int  setGeodParams ( lua_State *L );
int  printSolver   ( lua_State *L );

int  setLocalTetrad( lua_State *L );
int  printTetrad   ( lua_State *L );

int  calculateGeodesic ( lua_State *L );

int  calculateParTransport ( lua_State *L );

int  interpolatePosByCoord  ( lua_State *L );
int  interpolatePosByLambda ( lua_State *L );

int  getCurrentLambda  ( lua_State *L );
int  getCurrentPos     ( lua_State *L );
int  getCurrentDir     ( lua_State *L );

int  saveGeodesic      ( lua_State* L );

#endif // M4D_LUA_SOLVER_H
