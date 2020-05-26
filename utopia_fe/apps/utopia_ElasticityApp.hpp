#ifndef UTOPIA_ELASTICITY_APP
#define UTOPIA_ELASTICITY_APP

#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {

    class ElasticityApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "-elasticity"; }
    };
}  // namespace utopia

#endif  // UTOPIA_ELASTICITY_APP
