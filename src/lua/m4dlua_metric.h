/* -------------------------------------------------------------------------------
    m4dlua_metric.h

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

#ifndef M4D_LUA_METRIC_H
#define M4D_LUA_METRIC_H

#include "lua/m4dlua.h"
#include "lua/m4dlua_utils.h"

void lua_reg_metric(lua_State *L);

int  printMetricDB   ( lua_State* L );

int  setMetric       ( lua_State *L );
int  setMetricParam  ( lua_State *L );
int  setMetricParams ( lua_State *L );
int  getMetricParam  ( lua_State *L );
int  printMetric     ( lua_State *L );

int  localToCoord    ( lua_State *L );
int  coordToLocal    ( lua_State *L );

int  transToPseudoCoords ( lua_State *L );
int  coordTrans(lua_State *L);

int  calcTidalMatrix( lua_State *L );

#endif // M4D_LUA_METRIC_H
