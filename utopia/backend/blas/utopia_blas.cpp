#include "utopia_blas.hpp"
#include "utopia_BiCGStab_impl.hpp"
#include "utopia_blas_Vector.hpp"
#include "utopia_ConjugateGradient_impl.hpp"

//explicit instantiations
namespace utopia {
    template class BiCGStab<BlasMatrixd, BlasVectord, HOMEMADE>;
    template class ConjugateGradient<BlasMatrixd, BlasVectord>;
    template class ConjugateGradient<BlasMatrixd, BlasVectord, HOMEMADE>;
}
