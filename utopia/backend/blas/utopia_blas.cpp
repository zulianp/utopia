#include "utopia_blas.hpp"
#include "utopia_BiCGStab_impl.hpp"
#include "utopia_blas_Vector.hpp"

//explicit instantiations
namespace utopia {
    // template class Wrapper<utopia::BlasDenseMatrix<double>, 2>;
    // template class Wrapper<std::vector<double>, 1>;
    // template class Wrapper<utopia::CRSMatrix<double>, 2>;
    // template class Wrapper<utopia::CCSMatrix<double>, 2>;

    template class BiCGStab<utopia::BlasDenseMatrix<double>, Vectord, HOMEMADE>;
    // template class BiCGStab<utopia::CRSMatrix<double>, Vectord, HOMEMADE>;
    // template class BiCGStab<utopia::CCSMatrix<double>, Vectord, HOMEMADE>;
}
