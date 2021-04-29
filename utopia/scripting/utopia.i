/* utopia.i */
 %module utopia
 %{
#include "utopia_script.hpp"
#define SWIG_FILE_WITH_INIT
 %}

#ifdef UTOPIA_WITH_NUMPY

%include "numpy.i"

%init %{
    import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* seq, int n)};

#else

%include "carrays.i"
%array_functions(double, double_array);

%typemap(out) double *from_utopia_to_carray(int size) %{
  $result = PyList_New(9); 
  for (int i = 0; i < 9; ++i) {
    PyList_SetItem($result, i, PyFloat_FromDouble($1[i]));
  }
  delete $1;
%}

%init %{
%}

#endif //UTOPIA_WITH_NUMPY

#include "utopia_script.hpp"