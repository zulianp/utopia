/* utopia.i */
 %module utopia
 %{
#include "utopia_script.hpp"
#define SWIG_FILE_WITH_INIT
 %}

%include "carrays.i"
%array_functions(double, double_array);


%init %{
%}


#include "utopia_script.hpp"