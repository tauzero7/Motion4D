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

/**
 * @brief Register metric functions to be available from Lua
 * @param L  Lua state
 */
void lua_reg_metric(lua_State *L);

/**
 * @brief Print metric database
 * @param L  Lua state
 * @return
 */
int  printMetricDB   ( lua_State* L );

/**
 * @brief Set metric
 *
 *    SetMetric("Schwarzschild", {mass = 1.0})
 *
 * @param L  Lua state
 * @return
 */
int  setMetric( lua_State *L );

/**
 * @brief Set metric parameters
 * @param L  Lua state
 * @return
 */
int  setMetricParam  ( lua_State *L );
int  setMetricParams ( lua_State *L );

/**
 * @brief Get metric parameter
 * @param L  Lua state
 * @return
 */
int  getMetricParam  ( lua_State *L );

/**
 * @brief Print the currently set metric.
 * @param L
 * @return
 */
int printMetric( lua_State *L );

/**
 * @brief localToCoord transformation
 *
 *    cpos = LocalToCoord({
 *        pos = {t, x, y, z},
 *        vec = {dt, dx, dy, dz},
 *        type = "lnrf"
 *    })
 *
 * @param L
 * @return
 */
int  localToCoord    ( lua_State *L );

/**
 * @brief coordToLocal transformation
 *
 *    lpos = CoordToLocal({
 *        pos = {t, x, y, z},
 *        vec = {dt, dx, dy, dz},
 *        type = "lnrf"
 *    })
 *
 * @param L
 * @return
 */
int  coordToLocal    ( lua_State *L );

int  transToPseudoCoords ( lua_State *L );
int  coordTrans(lua_State *L);


/**
 * @brief Calculate product between two 4-vectors 'u' and 'v' given at a specific position 'pos'.
 *   Both vectors have to be given in coordinate representation.
 *
 *    prod = CalcProduct({
 *      pos = {t, x, y, z},
 *      u = {u0, u1, u2, u3},
 *      v = {v0, v1, v2, v3}
 *    })
 * @param L
 * @return
 */
int  calcProduct( lua_State *L );


int  calcTidalMatrix( lua_State *L );

#endif // M4D_LUA_METRIC_H
