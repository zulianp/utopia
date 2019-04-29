#ifndef UTOPIA_LEAST_SQUARES_HELMHOLTZ_HPP
#define UTOPIA_LEAST_SQUARES_HELMHOLTZ_HPP

#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {

    class LeastSquaresHelmholtzApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "-helmholtz";
        }
    };
}


#endif //UTOPIA_LEAST_SQUARES_HELMHOLTZ_HPP
