#ifndef UTOPIA_CONVERT_TO_DIEGO_MESH_APP_HPP
#define UTOPIA_CONVERT_TO_DIEGO_MESH_APP_HPP

#include "utopia_FEApp.hpp"
#include <string>


namespace utopia {

    class ConvertToDiegoMeshApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "-convert_to_diego";
        }
    };
}


#endif //UTOPIA_CONVERT_TO_DIEGO_MESH_APP_HPP
