#include "utopia_trilinos.hpp"

#include "utopia_BiCGStab_impl.hpp"
#include "utopia_ConjugateGradient_impl.hpp"
#include "utopia_Multigrid.hpp"
#include "utopia_Wrapper.hpp"

#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_Vector.hpp"
#include "utopia_trilinos_RowView.hpp"

#include "utopia_Tpetra_Matrix_impl.hpp"
#include "utopia_Tpetra_Vector_impl.hpp"

namespace utopia {

    template class RowView<TpetraMatrixd>;
    template class BiCGStab<TpetraMatrixd, TpetraVectord>;
    template class BiCGStab<TpetraMatrixd, TpetraVectord, HOMEMADE>;
    template class ConjugateGradient<TpetraMatrixd, TpetraVectord>;
    template class ConjugateGradient<TpetraMatrixd, TpetraVectord, HOMEMADE>;
    template class Factorization<TpetraMatrixd, TpetraVectord>;
    template class Multigrid<TpetraMatrixd, TpetraVectord, HOMEMADE>;

}  // namespace utopia
