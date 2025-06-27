// Mesh2pcfCpp.i

%module(directors="1") Mesh2pcfCpp

%include "stl.i"
%include "stdint.i"
%include "std_string.i"
%include "std_vector.i"
%include <std_shared_ptr.i>
%include <typemaps.i>
%include "std_iostream.i"
%include numpy.i

// Generic standard exceptions handling
%include "exception.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%{
#define SWIG_FILE_WITH_INIT

#include "Mesh2pcfCpp.h"
%}


%include "Mesh2pcfCpp.h"

%template(DoubleVector) std::vector<double>;
%template(DoubleVectorVector) std::vector<std::vector<double>>;
%template(DoubleGrid) std::vector<std::vector<std::vector<double>>>;
