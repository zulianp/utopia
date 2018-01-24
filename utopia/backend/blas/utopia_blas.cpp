#include "utopia_blas.hpp"

//explicit instantiations
namespace utopia {
	template class Wrapper<utopia::Matrix<double>, 2>;
	template class Wrapper<std::vector<double>, 1>;
	template class Wrapper<utopia::CRSMatrix<double>, 2>;
	template class Wrapper<utopia::CCSMatrix<double>, 2>;
}
