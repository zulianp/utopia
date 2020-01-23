#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
#include "utopia_trilinos_StructuredGrid.hpp"

namespace utopia {
    template class Mesh<UniformQuad4<double>, TrilinosCommunicator, TpetraVector::ExecutionSpace, Uniform<>>;
}

#endif //WITH_TRILINOS