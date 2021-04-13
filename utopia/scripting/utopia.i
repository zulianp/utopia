
/* utopia.i */
 %module utopia
 %{
#include "utopia_script.hpp"
#define SWIG_FILE_WITH_INIT
 %}

%include "numpy.i"


%init %{
    import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* seq, int n)};
#include "utopia_script.hpp"


