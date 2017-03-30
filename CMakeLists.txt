
cmake_minimum_required(VERSION 2.6)
project(libmotion4d)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR})

# set(GSL_DIR $ENV{HOME}/local/gsl/1.15)
set(GSL_DIR      $ENV{GSL_INSTALL_PATH}  CACHE FILEPATH "GSL_DIR")
set(GSL_LIB_DIR  $ENV{GSL_LIB_PATH}      CACHE FILEPATH "GSL_LIB_DIR")
#set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-rpath ${GSL_LIB_DIR}")

set(LUA_AVAILABLE OFF CACHE BOOL "have lua library available")
set(LUA_DIR       $ENV{LUA_HOME}        CACHE FILEPATH "Lua root path")
set(LUA_LIB_DIR   $ENV{LUA_HOME}/lib    CACHE FILEPATH "Lua library path")

set(COMPILE_TEST_PRGS OFF CACHE BOOL "Do compile all test programs")

set(M4D_ROOT_DIR   src)
set(M4D_EXTRA_DIR  src/extra)
set(M4D_MATH_DIR   src/math)
set(M4D_METRIC_DIR src/metric)
set(M4D_MOTION_DIR src/motion)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${GSL_DIR}/include)
link_directories(${GSL_LIB_DIR})

# Gather source files and store the result 
# in the variable 'm4d_source_files'.
file( GLOB m4d_source_files ${M4D_ROOT_DIR}/*.h 
    ${M4D_EXTRA_DIR}/*.cpp ${M4D_EXTRA_DIR}/*.h
    ${M4D_MATH_DIR}/*.cpp ${M4D_MATH_DIR}/*.h
    ${M4D_METRIC_DIR}/*.cpp ${M4D_METRIC_DIR}/*.h
    ${M4D_MOTION_DIR}/*.cpp ${M4D_MOTION_DIR}/*.h
    )
list( REMOVE_ITEM m4d_source_files 
    ${CMAKE_SOURCE_DIR}/${M4D_METRIC_DIR}/m4dMetricAlcubierreAccel.cpp )


if (LUA_AVAILABLE)
    set(M4D_LUA_DIR src/lua)
    include_directories(${LUA_DIR}/include)
    link_directories(${LUA_LIB_DIR})
    file( GLOB  lua_source_files ${M4D_LUA_DIR}/*.cpp ${M4D_LUA_DIR}/*.h)
    include_directories(${M4D_LUA_DIR})
endif()

# In debug mode an additional "d" is appended to any library
set(DAR "")
set(WDAR "")
if(CMAKE_BUILD_TYPE MATCHES "debug" OR 
    CMAKE_BUILD_TYPE MATCHES "Debug" OR 
    CMAKE_CONFIGURATION_TYPES MATCHES "debug" OR 
    CMAKE_CONFIGURATION_TYPES MATCHES "Debug")
    set(DAR "d")
    if(WIN32)
       set(WDAR "d")
    else()
      add_definitions(-pedantic)
      set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -DDEBUG -D_DEBUG -ggdb")
    endif(WIN32)
endif()
#message(m4d${DAR})

if(WIN32)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DM4D_LIB -D_m4d_EXPORTS")
endif(WIN32)

# build m4d library
add_library(m4d${DAR} SHARED ${m4d_source_files})
target_link_libraries(m4d${DAR} gsl gslcblas)
set_target_properties(m4d${DAR} 
    PROPERTIES  
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

install(TARGETS m4d${DAR} LIBRARY DESTINATION lib  RUNTIME DESTINATION bin)
install(TARGETS m4d${DAR} DESTINATION lib)
install(FILES src/m4dGlobalDefs.h DESTINATION include)
install(DIRECTORY ${M4D_EXTRA_DIR}  DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)
install(DIRECTORY ${M4D_MATH_DIR}   DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)
install(DIRECTORY ${M4D_METRIC_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE PATTERN kickedout EXCLUDE)
install(DIRECTORY ${M4D_MOTION_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)

# If also test programs shall be compiled
if (COMPILE_TEST_PRGS)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    #link_directories(${CMAKE_CACHEFILE_DIR}/Debug)
    
    # calculate geodesic
    set(TEST_CALC_GEOD_PATH  ${CMAKE_SOURCE_DIR}/test/calcGeodesic)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${TEST_CALC_GEOD_PATH})
    add_executable(m4dCalcGeodesic${DAR} ${TEST_CALC_GEOD_PATH}/calcGeodesic.cpp)
    target_link_libraries(m4dCalcGeodesic${DAR} m4d${DAR} gsl gslcblas)
    install(TARGETS m4dCalcGeodesic${DAR} DESTINATION bin)
    install(FILES ${TEST_CALC_GEOD_PATH}/kerr.ini DESTINATION bin)
    install(FILES ${TEST_CALC_GEOD_PATH}/plotGeodesic.gnu DESTINATION bin)

    # calculate parallel transport
    set(TEST_CALC_PAR_PATH ${CMAKE_SOURCE_DIR}/test/calcParallel)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${TEST_CALC_PAR_PATH})
    add_executable(m4dCalcParallel${DAR} ${TEST_CALC_PAR_PATH}/calcParallel.cpp)
    target_link_libraries(m4dCalcParallel${DAR} m4d${DAR} gsl gslcblas)
    install(TARGETS m4dCalcParallel${DAR} DESTINATION bin)
    install(FILES ${TEST_CALC_PAR_PATH}/kerr_timelike.ini DESTINATION bin)
    
    # calculate circular geodesic
    set(TEST_CIRC_GEOD_PATH ${CMAKE_SOURCE_DIR}/test/testCircular)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${TEST_CIRC_GEOD_PATH})
    add_executable(m4dTestCircular${DAR} ${TEST_CIRC_GEOD_PATH}/testCircular.cpp)
    target_link_libraries(m4dTestCircular${DAR} m4d${DAR} gsl gslcblas)
    install(TARGETS m4dTestCircular${DAR} DESTINATION bin)
    
    # print metric database
    set(TEST_METRIC_DATABASE_PATH ${CMAKE_SOURCE_DIR}/test/testDatabase)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${TEST_METRIC_DATABASE_PATH})
    add_executable(m4dTestDatabase${DAR} ${TEST_METRIC_DATABASE_PATH}/testDatabase.cpp)
    target_link_libraries(m4dTestDatabase${DAR} m4d${DAR} gsl gslcblas)
    install(TARGETS m4dTestDatabase${DAR} DESTINATION bin)

    # calculate fermi-walker transport
    set(TEST_FW_TRANSPORT_PATH ${CMAKE_SOURCE_DIR}/test/testFWworldline)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${TEST_FW_TRANSPORT_PATH})
    add_executable(m4dTestFWworldline${DAR} ${TEST_FW_TRANSPORT_PATH}/testFWworldline.cpp)
    target_link_libraries(m4dTestFWworldline${DAR} m4d${DAR} gsl gslcblas)
    install(TARGETS m4dTestFWworldline${DAR} DESTINATION bin)

    # test geodesic calculation
    set(TEST_GEOD_PATH ${CMAKE_SOURCE_DIR}/test/testGeodesic)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${TEST_GEOD_PATH})
    add_executable(m4dTestGeodesic${DAR} ${TEST_GEOD_PATH}/testGeodesic.cpp)
    target_link_libraries(m4dTestGeodesic${DAR} m4d${DAR} gsl gslcblas)
    install(TARGETS m4dTestGeodesic${DAR} DESTINATION bin)
    install(FILES ${TEST_GEOD_PATH}/plotSchwGeods.gnu DESTINATION bin)

    # test Jacobi calculation
    set(TEST_JACOBI_CALC_PATH ${CMAKE_SOURCE_DIR}/test/testJacobi)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${TEST_JACOBI_CALC_PATH})
    add_executable(m4dTestJacobi${DAR} ${TEST_JACOBI_CALC_PATH}/testJacobi.cpp)
    target_link_libraries(m4dTestJacobi${DAR} m4d${DAR} gsl gslcblas)
    install(TARGETS m4dTestJacobi${DAR} DESTINATION bin)
    install(FILES ${TEST_JACOBI_CALC_PATH}/kerr_lens.ini DESTINATION bin) 
    install(FILES ${TEST_JACOBI_CALC_PATH}/plotJacobi.gnu DESTINATION bin)  

    # test simple geodesic calculation
    set(TEST_SIMPLE_GEOD_PATH ${CMAKE_SOURCE_DIR}/test/testSimpleGeodesic)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${TEST_SIMPLE_GEOD_PATH})
    add_executable(m4dTestSimpleGeodesic${DAR} ${TEST_SIMPLE_GEOD_PATH}/testSimpleGeodesic.cpp)
    target_link_libraries(m4dTestSimpleGeodesic${DAR} m4d${DAR} gsl gslcblas)
    install(TARGETS m4dTestSimpleGeodesic${DAR} DESTINATION bin)
    
    # test cp motion calculation
    set(TEST_CP_MOTION_PATH ${CMAKE_SOURCE_DIR}/test/testCPmotion)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_CP_MOTION_PATH})
    add_executable(m4dTestCPmotion${DAR} ${TEST_CP_MOTION_PATH}/testCPmotion.cpp)
    target_link_libraries(m4dTestCPmotion${DAR} m4d${DAR} gsl gslcblas)
    
    install(FILES ${CMAKE_SOURCE_DIR}/test/testAll.bash DESTINATION bin)

    if(WIN32)
        add_custom_command(TARGET m4d${DAR} 
        POST_BUILD 
        COMMAND copy /Y $(LocalDebuggerWorkingDirectory)$(Configuration)\\$(TargetName).dll \"${TEST_CALC_GEOD_PATH}\"\\    
        COMMAND copy /Y $(LocalDebuggerWorkingDirectory)$(Configuration)\\$(TargetName).dll \"${TEST_CALC_PAR_PATH}\"\\    
        COMMAND copy /Y $(LocalDebuggerWorkingDirectory)$(Configuration)\\$(TargetName).dll \"${TEST_CIRC_GEOD_PATH}\"\\
        COMMAND copy /Y $(LocalDebuggerWorkingDirectory)$(Configuration)\\$(TargetName).dll \"${TEST_METRIC_DATABASE_PATH}\"\\
        COMMAND copy /Y $(LocalDebuggerWorkingDirectory)$(Configuration)\\$(TargetName).dll \"${TEST_FW_TRANSPORT_PATH}\"\\
        COMMAND copy /Y $(LocalDebuggerWorkingDirectory)$(Configuration)\\$(TargetName).dll \"${TEST_GEOD_PATH}\"\\
        COMMAND copy /Y $(LocalDebuggerWorkingDirectory)$(Configuration)\\$(TargetName).dll \"${TEST_JACOBI_CALC_PATH}\"\\
        COMMENT "Copy .dll into main folder")
    endif(WIN32)    

endif(COMPILE_TEST_PRGS)


if (LUA_AVAILABLE)
    # build m4d_lua_library
    add_library(m4d_lua${DAR} ${lua_source_files})
    target_link_libraries(m4d_lua${DAR} m4d${DAR} lua${WDAR})
    set_target_properties(m4d_lua${DAR} PROPERTIES  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
    install(TARGETS m4d_lua${DAR} DESTINATION lib)
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/lua DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)

    # build m4d lua program
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR})
    add_executable(m4dLua${DAR} lua/m4dlua_main.cpp ${lua_source_files})
    if(WIN32)
        target_link_libraries(m4dLua${DAR} m4d${DAR} lua${WDAR})
        add_custom_command(TARGET m4dLua${DAR} 
           POST_BUILD 
           COMMAND copy /Y $(LocalDebuggerWorkingDirectory)$(Configuration)\\m4d${DAR}.dll $(SolutionDir)..\\..\\$(Configuration)
        )
    else()
        target_link_libraries(m4dLua${DAR} m4d${DAR} m4d_lua${DAR} lua${WDAR} dl)
    endif(WIN32)
    install(TARGETS m4dLua${DAR} DESTINATION bin)
    install(FILES ${CMAKE_SOURCE_DIR}/lua/examples/schwarzschild.lua DESTINATION lua)  
endif()

