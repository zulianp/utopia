/* utopia.i */
 %module utopia
 %{
#include "utopia_script.hpp"
#define SWIG_FILE_WITH_INIT
 %}

#define UTOPIA_WITH_NUMPY
#ifdef UTOPIA_WITH_NUMPY

%include "numpy.i"

%init %{
    import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* seq, int n)};
%apply (double* INPLACE_ARRAY1, int DIM1) {(double*  seq_2, int n_2)};

#else

%include "carrays.i"
%array_functions(double, double_array);

%init %{
%}

#endif //UTOPIA_WITH_NUMPY

#include "utopia_script.hpp"
