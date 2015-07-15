
%module m4d

#pragma SWIG nowarn=362,389,503

%typemap(in) double*(double temp[4]) {   // temp[4] becomes a local variable
  int i;
  if (PyTuple_Check($input)) {
    if (!PyArg_ParseTuple($input,"dddd",temp,temp+1,temp+2,temp+3)) {
      PyErr_SetString(PyExc_TypeError,"tuple must have 4 elements");
      return NULL;
    }
    $1 = &temp[0];
  } else {
    PyErr_SetString(PyExc_TypeError,"expected a tuple.");
    return NULL;
  }
}

%typemap(out) std::tuple<double,double,double,double>  %{
    $result = PyTuple_New(4);

    PyTuple_SetItem($result, 0, PyFloat_FromDouble(std::get<0>($1)));
    PyTuple_SetItem($result, 1, PyFloat_FromDouble(std::get<1>($1)));
    PyTuple_SetItem($result, 2, PyFloat_FromDouble(std::get<2>($1)));
    PyTuple_SetItem($result, 3, PyFloat_FromDouble(std::get<3>($1)));
%}

%typemap(out) std::tuple<double,double,double>  %{
    $result = PyTuple_New(3);

    PyTuple_SetItem($result, 0, PyFloat_FromDouble(std::get<0>($1)));
    PyTuple_SetItem($result, 1, PyFloat_FromDouble(std::get<1>($1)));
    PyTuple_SetItem($result, 2, PyFloat_FromDouble(std::get<2>($1)));
%}


%{
#include "m4dGlobalDefs.h"
#include "VnD.h"
#include "Mat.h"
#include "m4dMetric.h"
#include "m4dMetricSchwarzschild.h"
#include "m4dMetricList.h"
#include "m4dMotion.h"
#include "m4dGeodesic.h"
#include "m4dGeodesicGSL.h"
#include "m4dGeodesicRK4.h"
#include "m4dObject.h"
%}


%include "../src/m4dGlobalDefs.h"
%include "../src/math/VnD.h"
%include "../src/math/Mat.h"
%include "../src/metric/m4dMetricList.h"
%include "../src/metric/m4dMetric.h"
%include "../src/metric/m4dMetricSchwarzschild.h"
%include "../src/metric/m4dMetricDatabase.h"
%include "../src/motion/m4dMotionList.h"
%include "../src/motion/m4dMotion.h"
%include "../src/motion/m4dGeodesic.h"
%include "../src/motion/m4dGeodesicGSL.h"
%include "../src/motion/m4dGeodesicRK4.h"
%include "../src/motion/m4dMotionDatabase.h"
%include "../src/extra/m4dObject.h"

%template(vec3) m4d::VnD<double,3>;
%template(vec4) m4d::VnD<double,4>;


%extend m4d::VnD<double,3> {
    double __getitem__(unsigned int i) {
        return (*($self))[i];
    }
}

%extend m4d::VnD<double,4> {
    double __getitem__(unsigned int i) {
        return (*($self))[i];
    }
}


%inline %{
std::tuple<double,double,double,double> getTuple(m4d::VnD<double,4> &v) {
    std::tuple<double, double, double, double> res;
    std::get<0>(res) = v.x(0);
    std::get<1>(res) = v.x(1);
    std::get<2>(res) = v.x(2);
    std::get<3>(res) = v.x(3);
    return res;
}

std::tuple<double, double, double> getThree() {
    std::tuple<double, double, double> res;
    std::get<0>(res) = 1.0; 
    std::get<1>(res) = 2.0;
    std::get<2>(res) = 3.0;
    return res;
}

std::vector<m4d::VnD<double,4>> getEmptyList() {
    return std::vector<m4d::VnD<double,4>>();
}

std::vector<double> getEmptyLambda() {
    return std::vector<double>();
}
%}

