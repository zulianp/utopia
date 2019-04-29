#ifndef UTOPIA_POISSON_APP_HPP
#define UTOPIA_POISSON_APP_HPP

#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {
    class PoissonApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "-poisson";
        }
    };
}


#endif //UTOPIA_POISSON_APP_HPP
