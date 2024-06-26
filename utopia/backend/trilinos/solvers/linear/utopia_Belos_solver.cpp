#include "utopia_Config.hpp"

#ifdef UTOPIA_ENABLE_TRILINOS_BELOS

#include "utopia_Belos_impl.hpp"
#include "utopia_trilinos_Types.hpp"

namespace utopia {
    template class BelosSolver<TpetraMatrixd, TpetraVectord>;
}  // namespace utopia

#endif  // UTOPIA_ENABLE_TRILINOS_BELOS
