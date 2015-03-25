
unix {
  GSL_DIR = # /usr/local/gsl/1.15
}

win32 {
  GSL_DIR = $$PWD/gsl
}

GSL_LIB_DIR = $$GSL_DIR/lib

CONFIG  += warn_on
QT -= gui core


LUA_DIR = /usr/local/lua/5.3
LUA_LIB_DIR = $$LUA_DIR/lib
CONFIG += USE_LUA


isEmpty(PREFIX) {
    PREFIX = /usr/local
}


p_global.path  = $$PREFIX/include
p_global.files = $$M4D_GLOBAL_HEADER
p_extra.path   = $$PREFIX/include/extra
p_extra.files  = $$M4D_EXTRA_HEADERS
p_math.path    = $$PREFIX/include/math
p_math.files   = $$M4D_MATH_HEADERS
p_metric.path  = $$PREFIX/include/metric
p_metric.files = $$M4D_METRIC_HEADERS
p_motion.path  = $$PREFIX/include/motion
p_motion.files = $$M4D_MOTION_HEADERS
p_lib.path     = $$PREFIX/lib
p_lib.files    = lib/lib$$TARGET$$DEBUG.so.* lib/lib$$TARGET$$DEBUG.a
p_other.path   = $$PREFIX/bin
p_other.files  = test/testAll.bash

INSTALLS += p_global p_extra p_math p_metric p_motion p_lib p_other
