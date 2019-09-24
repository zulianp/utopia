#include "utopia_blas.hpp"
#include "utopia_BiCGStab_impl.hpp"
#include "utopia_blas_Vector.hpp"

//explicit instantiations
namespace utopia {
    template class BiCGStab<BlasMatrixd, BlasVectord, HOMEMADE>;
}
