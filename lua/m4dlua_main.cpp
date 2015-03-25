/* -------------------------------------------------------------------------------
    m4dlua_main.cpp

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

#include <iostream>

#include "lua/m4dlua.h"
#include "lua/m4dlua_metric.h"
#include "lua/m4dlua_solver.h"
#include "lua/m4dlua_utils.h"

lua_State*  L = NULL;

m4d::Object   mObject;



int main(int argc, char* argv[]) {
    if (argc!=2) {
        fprintf(stderr,"Usage: ./m4dLua <file.lua>\n");
        return -1;
    }
 
    L = luaL_newstate();
    luaL_openlibs(L);
    
    // register all functions
    lua_reg_metric(L);
    lua_reg_solver(L);
    lua_reg_utils(L);


    //int numErrors = 0;
    //if (luaL_loadfile(L, argv[1]) || lua_pcall(L, 0, 0, 0)) {
    int errornr = LUA_OK;
    if ((errornr = luaL_dofile(L, argv[1]))) {
        mlua_error(L, "Error: %s\n", lua_tostring(L, -1));
    } else {
        lua_close(L);
    }
    return 0;   
}
