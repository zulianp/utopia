#include "utopia_trilinos.hpp"
#include "utopia_BiCGStab_impl.hpp"
#include "utopia_ConjugateGradient_impl.hpp"
#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_Vector.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_trilinos_RowView.hpp"

namespace utopia {

    template class RowView<TpetraMatrixd>;
    template class BiCGStab<TpetraMatrixd, TpetraVectord, HOMEMADE>;
    template class ConjugateGradient<TpetraMatrixd, TpetraVectord>;
    template class ConjugateGradient<TpetraMatrixd, TpetraVectord, HOMEMADE>;

}  // namespace utopia
