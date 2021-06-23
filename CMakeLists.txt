# -----------------------------------------------------------------------------
#   CMakeLists.txt  for libMotion4D
#
#   Run cmake-gui to select a generator.
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.15)
project(libmotion4d)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR})

set(COMPILE_ALL_METRICS  OFF  CACHE  BOOL  "Compile all metrics")


# ---------------------------------------------
#  GNU Scientific Library (GSL)
# ---------------------------------------------
set(GSL_DIR      $ENV{GSL_INSTALL_PATH}  CACHE FILEPATH "GSL_DIR")
set(GSL_LIB_DIR  $ENV{GSL_LIB_PATH}      CACHE FILEPATH "GSL_LIB_DIR")

# ---------------------------------------------
#  LUA
# ---------------------------------------------
set(LUA_AVAILABLE OFF CACHE BOOL "have lua library available")
set(LUA_DIR       $ENV{LUA_HOME}        CACHE FILEPATH "Lua root path")
set(LUA_LIB_DIR   $ENV{LUA_HOME}/lib    CACHE FILEPATH "Lua library path")

# ---------------------------------------------
set(COMPILE_TEST_PRGS OFF CACHE BOOL "Do compile all test programs")

# ---------------------------------------------
#  source dirs
# ---------------------------------------------
set(M4D_ROOT_DIR   src)
set(M4D_EXTRA_DIR  src/extra)
set(M4D_MATH_DIR   src/math)
set(M4D_METRIC_DIR src/metric)
set(M4D_MOTION_DIR src/motion)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${GSL_DIR}/include)
link_directories(${GSL_LIB_DIR} lib)

# ---------------------------------------------
# Gather source files and store the result 
# in the variable 'm4d_source_files'.
# ---------------------------------------------
file( GLOB m4d_source_files ${M4D_ROOT_DIR}/*.h 
    ${M4D_EXTRA_DIR}/*.cpp ${M4D_EXTRA_DIR}/*.h
    ${M4D_MATH_DIR}/*.cpp ${M4D_MATH_DIR}/*.h
    ${M4D_METRIC_DIR}/*.cpp ${M4D_METRIC_DIR}/*.h
    ${M4D_MOTION_DIR}/*.cpp ${M4D_MOTION_DIR}/*.h
    )
list( REMOVE_ITEM m4d_source_files 
    ${CMAKE_SOURCE_DIR}/${M4D_METRIC_DIR}/m4dMetricAlcubierreAccel.cpp )


if (COMPILE_ALL_METRICS) 
    add_definitions(-DALL_METRICS_AVAILABLE)
endif()

if (LUA_AVAILABLE)
    set(M4D_LUA_DIR src/lua)
    include_directories(${LUA_DIR}/include)
    link_directories(${LUA_LIB_DIR})
    file( GLOB  lua_source_files ${M4D_LUA_DIR}/*.cpp ${M4D_LUA_DIR}/*.h)
    include_directories(${M4D_LUA_DIR})
endif()

# be pedantic
#add_definitions(-pedantic)


if(WIN32)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DM4D_LIB -D_m4d_EXPORTS")
endif(WIN32)

# ---------------------------------------------
#   build m4d library
# ---------------------------------------------
add_library(m4d SHARED ${m4d_source_files})
target_link_libraries(m4d gsl gslcblas)
set_target_properties(m4d PROPERTIES  
    DEBUG_POSTFIX "_d"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_SOURCE_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_SOURCE_DIR}/lib")

install(TARGETS m4d LIBRARY DESTINATION lib  RUNTIME DESTINATION bin)
install(TARGETS m4d DESTINATION lib)
install(FILES src/m4dGlobalDefs.h DESTINATION include)
install(DIRECTORY ${M4D_EXTRA_DIR}  DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)
install(DIRECTORY ${M4D_MATH_DIR}   DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)
install(DIRECTORY ${M4D_METRIC_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE PATTERN kickedout EXCLUDE)
install(DIRECTORY ${M4D_MOTION_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)

# ---------------------------------------------
#  If also test programs shall be compiled
# ---------------------------------------------
if (COMPILE_TEST_PRGS)
    #set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    
    # calculate geodesic
    set(TEST_CALC_GEOD_PATH  ${CMAKE_SOURCE_DIR}/test/calcGeodesic)    
    add_executable(m4dCalcGeodesic ${TEST_CALC_GEOD_PATH}/calcGeodesic.cpp)
    add_dependencies(m4dCalcGeodesic m4d)
    set_target_properties(m4dCalcGeodesic PROPERTIES  
        DEBUG_POSTFIX "_d"
        RUNTIME_OUTPUT_DIRECTORY  ${TEST_CALC_GEOD_PATH}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${TEST_CALC_GEOD_PATH}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${TEST_CALC_GEOD_PATH})
    target_link_libraries(m4dCalcGeodesic m4d gsl gslcblas)
    install(TARGETS m4dCalcGeodesic DESTINATION bin)
    install(FILES ${TEST_CALC_GEOD_PATH}/kerr.ini DESTINATION bin)
    install(FILES ${TEST_CALC_GEOD_PATH}/plotGeodesic.gnu DESTINATION bin)

    # calculate parallel transport
    set(TEST_CALC_PAR_PATH ${CMAKE_SOURCE_DIR}/test/calcParallel)
    add_executable(m4dCalcParallel ${TEST_CALC_PAR_PATH}/calcParallel.cpp)
    add_dependencies(m4dCalcParallel m4d)
    set_target_properties(m4dCalcParallel PROPERTIES  
        DEBUG_POSTFIX "_d"
        RUNTIME_OUTPUT_DIRECTORY  ${TEST_CALC_PAR_PATH}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${TEST_CALC_PAR_PATH}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${TEST_CALC_PAR_PATH})
    target_link_libraries(m4dCalcParallel m4d gsl gslcblas)
    install(TARGETS m4dCalcParallel DESTINATION bin)
    install(FILES ${TEST_CALC_PAR_PATH}/kerr_timelike.ini DESTINATION bin)
    
    # calculate circular geodesic
    set(TEST_CIRC_GEOD_PATH ${CMAKE_SOURCE_DIR}/test/testCircular)
    add_executable(m4dTestCircular ${TEST_CIRC_GEOD_PATH}/testCircular.cpp)
    add_dependencies(m4dTestCircular m4d)
    set_target_properties(m4dTestCircular PROPERTIES  
        DEBUG_POSTFIX "_d"
        RUNTIME_OUTPUT_DIRECTORY  ${TEST_CIRC_GEOD_PATH}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${TEST_CIRC_GEOD_PATH}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${TEST_CIRC_GEOD_PATH})
    target_link_libraries(m4dTestCircular m4d gsl gslcblas)
    install(TARGETS m4dTestCircular DESTINATION bin)
    
    # print metric database
    set(TEST_METRIC_DATABASE_PATH ${CMAKE_SOURCE_DIR}/test/testDatabase)
    add_executable(m4dTestDatabase ${TEST_METRIC_DATABASE_PATH}/testDatabase.cpp)
    add_dependencies(m4dTestDatabase m4d)
    set_target_properties(m4dTestDatabase PROPERTIES  
        DEBUG_POSTFIX "_d"
        RUNTIME_OUTPUT_DIRECTORY  ${TEST_METRIC_DATABASE_PATH}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${TEST_METRIC_DATABASE_PATH}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${TEST_METRIC_DATABASE_PATH})
    target_link_libraries(m4dTestDatabase m4d gsl gslcblas)
    install(TARGETS m4dTestDatabase DESTINATION bin)

    # calculate fermi-walker transport
    set(TEST_FW_TRANSPORT_PATH ${CMAKE_SOURCE_DIR}/test/testFWworldline)
    add_executable(m4dTestFWworldline ${TEST_FW_TRANSPORT_PATH}/testFWworldline.cpp)
    add_dependencies(m4dTestFWworldline m4d)
    set_target_properties(m4dTestFWworldline PROPERTIES  
        DEBUG_POSTFIX "_d"
        RUNTIME_OUTPUT_DIRECTORY  ${TEST_FW_TRANSPORT_PATH}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${TEST_FW_TRANSPORT_PATH}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${TEST_FW_TRANSPORT_PATH})
    target_link_libraries(m4dTestFWworldline m4d gsl gslcblas)
    install(TARGETS m4dTestFWworldline DESTINATION bin)

    # test geodesic calculation
    set(TEST_GEOD_PATH ${CMAKE_SOURCE_DIR}/test/testGeodesic)
    add_executable(m4dTestGeodesic ${TEST_GEOD_PATH}/testGeodesic.cpp)
    add_dependencies(m4dTestGeodesic m4d)
    set_target_properties(m4dTestGeodesic PROPERTIES  
        DEBUG_POSTFIX "_d"
        RUNTIME_OUTPUT_DIRECTORY  ${TEST_GEOD_PATH}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${TEST_GEOD_PATH}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${TEST_GEOD_PATH})
    target_link_libraries(m4dTestGeodesic m4d gsl gslcblas)
    install(TARGETS m4dTestGeodesic DESTINATION bin)
    install(FILES ${TEST_GEOD_PATH}/plotSchwGeods.gnu DESTINATION bin)

    # test Jacobi calculation
    set(TEST_JACOBI_CALC_PATH ${CMAKE_SOURCE_DIR}/test/testJacobi)
    add_executable(m4dTestJacobi ${TEST_JACOBI_CALC_PATH}/testJacobi.cpp)
    add_dependencies(m4dTestJacobi m4d)
    set_target_properties(m4dTestJacobi PROPERTIES  
        DEBUG_POSTFIX "_d"
        RUNTIME_OUTPUT_DIRECTORY  ${TEST_JACOBI_CALC_PATH}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${TEST_JACOBI_CALC_PATH}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${TEST_JACOBI_CALC_PATH})
    target_link_libraries(m4dTestJacobi m4d gsl gslcblas)
    install(TARGETS m4dTestJacobi DESTINATION bin)
    install(FILES ${TEST_JACOBI_CALC_PATH}/kerr_lens.ini DESTINATION bin) 
    install(FILES ${TEST_JACOBI_CALC_PATH}/plotJacobi.gnu DESTINATION bin)  

    # test simple geodesic calculation
    set(TEST_SIMPLE_GEOD_PATH ${CMAKE_SOURCE_DIR}/test/testSimpleGeodesic)
    add_executable(m4dTestSimpleGeodesic ${TEST_SIMPLE_GEOD_PATH}/testSimpleGeodesic.cpp)
    add_dependencies(m4dTestSimpleGeodesic m4d)
    set_target_properties(m4dTestSimpleGeodesic PROPERTIES  
        DEBUG_POSTFIX "_d"
        RUNTIME_OUTPUT_DIRECTORY  ${TEST_SIMPLE_GEOD_PATH}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${TEST_SIMPLE_GEOD_PATH}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${TEST_SIMPLE_GEOD_PATH})
    target_link_libraries(m4dTestSimpleGeodesic m4d gsl gslcblas)
    install(TARGETS m4dTestSimpleGeodesic DESTINATION bin)
    
    # test cp motion calculation
    set(TEST_CP_MOTION_PATH ${CMAKE_SOURCE_DIR}/test/testCPmotion)
    add_executable(m4dTestCPmotion ${TEST_CP_MOTION_PATH}/testCPmotion.cpp)
    add_dependencies(m4dTestCPmotion m4d)
    set_target_properties(m4dTestCPmotion PROPERTIES  
        DEBUG_POSTFIX "_d"
        RUNTIME_OUTPUT_DIRECTORY  ${TEST_CP_MOTION_PATH}
        RUNTIME_OUTPUT_DIRECTORY_DEBUG  ${TEST_CP_MOTION_PATH}
        RUNTIME_OUTPUT_DIRECTORY_RELEASE  ${TEST_CP_MOTION_PATH})
    target_link_libraries(m4dTestCPmotion m4d gsl gslcblas)
    
    install(FILES ${CMAKE_SOURCE_DIR}/test/testAll.bash DESTINATION bin)

    if(WIN32)
        add_custom_command(TARGET m4d
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

# ---------------------------------------------
#  build m4d_lua_library
# ---------------------------------------------
if (LUA_AVAILABLE)
    add_library(m4d_lua${DAR} ${lua_source_files})
    target_link_libraries(m4d_lua${DAR} m4d lua${WDAR})
    set_target_properties(m4d_lua${DAR} PROPERTIES  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
    install(TARGETS m4d_lua${DAR} DESTINATION lib)
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/lua DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)

    # build m4d lua program
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR})
    add_executable(m4dLua${DAR} lua/m4dlua_main.cpp ${lua_source_files})
    if(WIN32)
        target_link_libraries(m4dLua${DAR} m4d lua${WDAR})
        add_custom_command(TARGET m4dLua${DAR} 
           POST_BUILD 
           COMMAND copy /Y $(LocalDebuggerWorkingDirectory)$(Configuration)\\m4d${DAR}.dll $(SolutionDir)..\\..\\$(Configuration)
        )
    else()
        target_link_libraries(m4dLua${DAR} m4d m4d_lua${DAR} lua${WDAR} dl)
    endif(WIN32)
    install(TARGETS m4dLua${DAR} DESTINATION bin)
    install(FILES ${CMAKE_SOURCE_DIR}/lua/examples/schwarzschild.lua DESTINATION lua)  
endif()

