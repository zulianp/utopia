#ifndef UTOPIA_FE_BACKEND_HPP
#define UTOPIA_FE_BACKEND_HPP

#include "utopia_Traits.hpp"

namespace utopia {

    template<int BackendType>
    class FEBackend {
    public:
        FEBackend()
        {
            static_assert(BackendType < HOMEMADE, "No Backend implemented");
        }
    };
}

#endif //UTOPIA_FE_BACKEND_HPP
