
include(../../common.pri)


M4D_DIR     = $$PWD/../../
TEMPLATE = app

CONFIG(debug, debug|release) {
    TARGET      = m4dCalcParallel_debug
    DESTDIR     = .
    OBJECTS_DIR = $$M4D_DIR/compiled/debug
    unix {
       LIBS += -Wl,-rpath $$M4D_DIR/lib -L$$M4D_DIR/lib -lm4d_debug
    }
}

CONFIG(release, debug|release) {
    TARGET      = m4dCalcParallel
    DESTDIR     = .
    OBJECTS_DIR = $$M4D_DIR/compiled/release
    unix {
       LIBS += -Wl,-rpath $$M4D_DIR/lib -L$$M4D_DIR/lib -lm4d
    }
}

INCLUDEPATH += . $$M4D_DIR $$M4D_DIR/src

SOURCES = calcParallel.cpp

unix {
    LIBS += -Wl,-rpath $$GSL_LIB_DIR -L$$GSL_LIB_DIR -lgsl -lgslcblas
    USE_LUA {
        LIBS += -L$$LUA_LIB_DIR -llua -ldl
    }
}

p_app.path  = $$PREFIX/bin
p_app.files = $$TARGET
p_other.path  = $$PREFIX/bin
p_other.files = kerr_timelike.ini plotGeodesics.gnu

INSTALLS += p_app p_other
