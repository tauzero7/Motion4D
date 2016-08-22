
######################################################################
M4D_DIR     = $$PWD
M4D_SRC_DIR = $$M4D_DIR/src

include(m4d_sources.pri)
include(m4d_lua.pri)
include(common.pri)

######################################################################
CONFIG  += warn_on
QT -= gui core
TEMPLATE = lib
##CONFIG += staticlib

unix {
QMAKE_CXXFLAGS += -std=c++11
}

CONFIG(debug, debug|release) {
    TARGET = m4d_debug
    DESTDIR     = $$M4D_DIR/lib
    OBJECTS_DIR = $$M4D_DIR/compiled/debug
}

CONFIG(release, debug|release) {
    TARGET = m4d
    DESTDIR     = $$M4D_DIR/lib
    OBJECTS_DIR = $$M4D_DIR/compiled/release
}


######################################################################  
HEADERS +=  $$M4D_HEADERS
SOURCES +=  $$M4D_SOURCES
INCLUDEPATH +=  . .. $$M4D_DIR $$M4D_SRC_DIR $$GSL_DIR/include

USE_LUA {
   HEADERS += $$LUA_HEADERS
   SOURCES += $$LUA_SOURCES
   INCLUDEPATH += $$LUA_DIR/include
   #DEFINES +=
}

unix {
    LIBS +=  -Wl,-rpath $$GSL_LIB_DIR\
              -L$$GSL_LIB_DIR -lgsl -lgslcblas
}

win32 {
    DEFINES += m4d_EXPORTS
    LIBS += -L"$$GSL_LIB_DIR"  $$GSL_LIB_DIR/gsl.lib  $$GSL_LIB_DIR/gslcblas.lib
    INCLUDEPATH += .  .. "$$M4D_DIR" "$$GSL_DIR/include"
}

