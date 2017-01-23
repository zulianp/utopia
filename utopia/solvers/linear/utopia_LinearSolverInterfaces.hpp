#ifndef UTOPIA_LINEAR_SOLVERS_INTERFACES_HPP
#define UTOPIA_LINEAR_SOLVERS_INTERFACES_HPP 

#include "utopia_ForwardDeclarations.hpp"

namespace utopia {
	template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
	class BiCGStab {};

	template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
	class GMRES {};

	template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
	class MINRES {};

	template<typename Matrix, typename Vector, int Backend = Traits<Vector>::Backend>
	class Factorization;

	template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
	class LUDecomposition {};

	template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
	class CholeskyDecomposition {};
}

#endif //UTOPIA_LINEAR_SOLVERS_INTERFACES_HPP
