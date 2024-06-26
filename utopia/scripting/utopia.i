/* utopia.i */
 %module utopia
 %{
#include "utopia_script.hpp"
#define SWIG_FILE_WITH_INIT
 %}

#define UTOPIA_ENABLE_NUMPY
#ifdef UTOPIA_ENABLE_NUMPY

%include "numpy.i"

%init %{
    import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* seq, int n)};

#else

%include "carrays.i"
%array_functions(double, double_array);

%init %{
%}

#endif //UTOPIA_ENABLE_NUMPY

#include "utopia_script.hpp"
