
#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS_AMESOS2

#include "utopia_Amesos2_impl.hpp"
#include "utopia_trilinos_Types.hpp"

namespace utopia {
    template class Amesos2Solver<TpetraMatrixd, TpetraVectord>;
}  // namespace utopia

#endif  // WITH_TRILINOS_AMESOS2