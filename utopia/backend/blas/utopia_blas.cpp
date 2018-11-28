#include "utopia_blas.hpp"
#include "utopia_BiCGStab_impl.hpp"

//explicit instantiations
namespace utopia {
	template class Wrapper<utopia::Matrix<double>, 2>;
	template class Wrapper<std::vector<double>, 1>;
	template class Wrapper<utopia::CRSMatrix<double>, 2>;
	template class Wrapper<utopia::CCSMatrix<double>, 2>;

	// template class BiCGStab<utopia::Matrix<double>, Vectord, HOMEMADE>;
	// template class BiCGStab<utopia::CRSMatrix<double>, Vectord, HOMEMADE>;
	// template class BiCGStab<utopia::CCSMatrix<double>, Vectord, HOMEMADE>;
}
