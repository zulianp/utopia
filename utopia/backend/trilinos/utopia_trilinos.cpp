#include "utopia_trilinos.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_Vector.hpp"
#include "utopia_trilinos_RowView.hpp"
#include "utopia_BiCGStab_impl.hpp"

namespace utopia{

    template class RowView<TpetraMatrix>;
    template class BiCGStab<TpetraMatrix, TpetraVector, HOMEMADE>;

}
