
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

set(COMPILE_TEST_PRGS OFF CACHE BOOL "do compile all test programs")

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
file( GLOB m4d_source_files ${M4D_ROOT_DIR}/*.cpp ${M4D_EXTRA_DIR}/*.cpp ${M4D_MATH_DIR}/*.cpp ${M4D_METRIC_DIR}/*.cpp ${M4D_MOTION_DIR}/*.cpp )
list( REMOVE_ITEM m4d_source_files ${CMAKE_SOURCE_DIR}/${M4D_METRIC_DIR}/m4dMetricAlcubierreAccel.cpp )


if (LUA_AVAILABLE)
    set(M4D_LUA_DIR src/lua)
    include_directories(${LUA_DIR}/include)
    link_directories(${LUA_LIB_DIR})
    file( GLOB  lua_source_files ${M4D_LUA_DIR}/*.cpp )
    include_directories(${M4D_LUA_DIR})
endif()

#set(BITS 32)
#if(CMAKE_SIZEOF_VOID_P EQUAL 8)
#  set(BITS 64)
#endif()

set(DAR "")
if(CMAKE_BUILD_TYPE MATCHES "debug" OR CMAKE_BUILD_TYPE MATCHES "DEBUG")
    set(DAR "_debug")
    add_definitions(-pedantic)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -DDEBUG -D_DEBUG -ggdb")
endif()

add_library(m4d${BITS}${DAR} SHARED ${m4d_source_files})
target_link_libraries(m4d${BITS}${DAR} gsl gslcblas)
set_target_properties(m4d${BITS}${DAR} PROPERTIES  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

install(TARGETS m4d${BITS}${DAR} DESTINATION lib)
install(FILES src/m4dGlobalDefs.h DESTINATION include)
install(DIRECTORY ${M4D_EXTRA_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)
install(DIRECTORY ${M4D_MATH_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)
install(DIRECTORY ${M4D_METRIC_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE PATTERN kickedout EXCLUDE)
install(DIRECTORY ${M4D_MOTION_DIR} DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)


if (COMPILE_TEST_PRGS)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR}/test/calcGeodesic)
    add_executable(m4dCalcGeodesic${DAR} test/calcGeodesic/calcGeodesic.cpp)
    target_link_libraries(m4dCalcGeodesic${DAR} m4d${BITS}${DAR} gsl gslcblas)
    install(TARGETS m4dCalcGeodesic${DAR} DESTINATION bin)
    install(FILES ${CMAKE_SOURCE_DIR}/test/calcGeodesic/kerr.ini DESTINATION bin)
    install(FILES ${CMAKE_SOURCE_DIR}/test/calcGeodesic/plotGeodesic.gnu DESTINATION bin)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR}/test/calcParallel)
    add_executable(m4dCalcParallel${DAR} test/calcParallel/calcParallel.cpp)
    target_link_libraries(m4dCalcParallel${DAR} m4d${BITS}${DAR} gsl gslcblas)
    install(TARGETS m4dCalcParallel${DAR} DESTINATION bin)
    install(FILES ${CMAKE_SOURCE_DIR}/test/calcParallel/kerr_timelike.ini DESTINATION bin)
    
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR}/test/testCircular)
    add_executable(m4dTestCircular${DAR} test/testCircular/testCircular.cpp)
    target_link_libraries(m4dTestCircular${DAR} m4d${BITS}${DAR} gsl gslcblas)
    install(TARGETS m4dTestCircular${DAR} DESTINATION bin)
    
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR}/test/testDatabase)
    add_executable(m4dTestDatabase${DAR} test/testDatabase/testDatabase.cpp)
    target_link_libraries(m4dTestDatabase${DAR} m4d${BITS}${DAR} gsl gslcblas)
    install(TARGETS m4dTestDatabase${DAR} DESTINATION bin)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR}/test/testFWworldline)
    add_executable(m4dTestFWworldline${DAR} test/testFWworldline/testFWworldline.cpp)
    target_link_libraries(m4dTestFWworldline${DAR} m4d${BITS}${DAR} gsl gslcblas)
    install(TARGETS m4dTestFWworldline${DAR} DESTINATION bin)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR}/test/testGeodesic)
    add_executable(m4dTestGeodesic${DAR} test/testGeodesic/testGeodesic.cpp)
    target_link_libraries(m4dTestGeodesic${DAR} m4d${BITS}${DAR} gsl gslcblas)
    install(TARGETS m4dTestGeodesic${DAR} DESTINATION bin)
    install(FILES ${CMAKE_SOURCE_DIR}/test/testGeodesic/plotSchwGeods.gnu DESTINATION bin)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR}/test/testSimpleGeodesic)
    add_executable(m4dTestSimpleGeodesic${DAR} test/testSimpleGeodesic/testSimpleGeodesic.cpp)
    target_link_libraries(m4dTestSimpleGeodesic${DAR} m4d${BITS}${DAR} gsl gslcblas)
    install(TARGETS m4dTestSimpleGeodesic${DAR} DESTINATION bin)

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR}/test/testJacobi)
    add_executable(m4dTestJacobi${DAR} test/testJacobi/testJacobi.cpp)
    target_link_libraries(m4dTestJacobi${DAR} m4d${BITS}${DAR} gsl gslcblas)
    install(TARGETS m4dTestJacobi${DAR} DESTINATION bin)
    install(FILES ${CMAKE_SOURCE_DIR}/test/testJacobi/kerr_lens.ini DESTINATION bin) 
    install(FILES ${CMAKE_SOURCE_DIR}/test/testJacobi/plotJacobi.gnu DESTINATION bin)  

    install(FILES ${CMAKE_SOURCE_DIR}/test/testAll.bash DESTINATION bin)
endif()


if (LUA_AVAILABLE)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_SOURCE_DIR})
	add_executable(m4dLua${DAR} lua/m4dlua_main.cpp ${lua_source_files})
	target_link_libraries(m4dLua${DAR} m4d${BITS}${DAR} lua dl)
	install(TARGETS m4dLua${DAR} DESTINATION bin)

	add_library(m4d_lua${BITS}${DAR} SHARED ${lua_source_files})
	set_target_properties(m4d_lua${BITS}${DAR} PROPERTIES  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
	install(TARGETS m4d_lua${BITS}${DAR} DESTINATION lib)
	install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/lua DESTINATION include FILES_MATCHING PATTERN "*.h" PATTERN .svn EXCLUDE PATTERN .deps EXCLUDE)
endif()

