 /* utopia.i */
 %module utopia
 %{
#include "utopia_script.hpp"
 %}

#include "utopia_script.hpp"

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}

%apply(double* IN_ARRAY1[ANY]) {(double * values)};