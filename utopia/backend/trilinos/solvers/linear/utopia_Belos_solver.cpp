#include <utopia_Config.hpp>

#ifdef WITH_TRILINOS_BELOS

#include "utopia_Belos_impl.hpp"
#include "utopia_trilinos_Types.hpp"

namespace utopia {
    template class BelosSolver<TpetraMatrixd, TpetraVectord>;
}
#endif  // WITH_TRILINOS_BELOS
