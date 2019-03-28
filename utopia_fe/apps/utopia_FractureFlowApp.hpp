#ifndef UTOPIA_FRACTURE_FLOW_APP_HPP
#define UTOPIA_FRACTURE_FLOW_APP_HPP

#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {
    class FractureFlowApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "-fracflow";
        }
    };
}

#endif //UTOPIA_FRACTURE_FLOW_APP_HPP
