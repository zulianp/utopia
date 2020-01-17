 /* utopia.i */
 %module utopia
 %{
#include "utopia_script.hpp"
 %}

 namespace algebra {
    class SparseMatrix {
    public:
        SparseMatrix();
        ~SparseMatrix();
    };
}
