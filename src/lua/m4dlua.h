/* -------------------------------------------------------------------------------
    m4dlua.h

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

#ifndef M4D_LUA_H
#define M4D_LUA_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

extern "C" {
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}

#include <m4dGlobalDefs.h>
#include <metric/m4dMetricDatabase.h>
#include <extra/m4dObject.h>
#include <extra/m4dUtilities.h>

extern m4d::Object  mObject;


#ifndef DEF_LUA_SIGNUM
#define DEF_LUA_SIGNUM(x) (x>=0 ? 1.0 : -1.0)
#endif

#ifndef DEF_LUA_MIN
#define DEF_LUA_MIN(x,y) ((x)<(y) ? (x) : (y))
#endif

#endif
