#ifndef UTOPIA_CONVERT_TO_DIEGO_MESH_APP_HPP
#define UTOPIA_CONVERT_TO_DIEGO_MESH_APP_HPP

#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {

    class ConvertToDiegoMeshApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "-convert_to_diego"; }
    };
}  // namespace utopia

#endif  // UTOPIA_CONVERT_TO_DIEGO_MESH_APP_HPP
