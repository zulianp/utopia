#ifndef UTOPIA_POLYMORPHIC_QP_SOLVER_HPP
#define UTOPIA_POLYMORPHIC_QP_SOLVER_HPP

#include "utopia_QPSolver.hpp"

namespace utopia {

	template<class Matrix, class Vector>
	class PolymorphicQPSolver : public QPSolver<Matrix, Vector>  {

	};

}

#endif //UTOPIA_POLYMORPHIC_QP_SOLVER_HPP
