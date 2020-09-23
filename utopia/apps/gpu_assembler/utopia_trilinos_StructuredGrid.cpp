#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_TRILINOS
#include "utopia_trilinos_StructuredGrid.hpp"

namespace utopia {
    template class Mesh<UniformQuad4<double>, TrilinosCommunicator, TpetraVector::ExecutionSpace, Uniform<>>;
}  // namespace utopia

#endif  // UTOPIA_WITH_TRILINOS