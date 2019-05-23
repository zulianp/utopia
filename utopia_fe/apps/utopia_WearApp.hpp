#ifndef UTOPIA_WEAR_APP
#define UTOPIA_WEAR_APP

#include "utopia_FEApp.hpp"
#include <string>

namespace utopia {

    class WearApp final : public FEApp {
    public:

        class SimulationInput;
        class SimulationOutput;

        void run(Input &in) override;

        inline static std::string command()
        {
            return "-wear";
        }
    };
}

#endif //UTOPIA_WEAR_APP
