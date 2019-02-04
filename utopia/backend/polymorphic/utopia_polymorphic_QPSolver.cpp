

#include "utopia_polymorphic_QPSolver_impl.hpp"
#include "utopia_Base.hpp"
#include "utopia_petsc_Types.hpp"


namespace utopia {

#ifdef WITH_PETSC
		template class PolymorphicQPSolver<DSMatrixd, DVectord>;
#endif //WITH_PETC

}