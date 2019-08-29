#ifndef UTOPIA_MIXED_POISSON_APP_HPP
#define UTOPIA_MIXED_POISSON_APP_HPP

#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {
    class MixedPoissonApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "-mixed_poisson";
        }
    };
}


#endif //UTOPIA_MIXED_POISSON_APP_HPP
