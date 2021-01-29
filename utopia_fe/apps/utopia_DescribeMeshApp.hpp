#ifndef UTOPIA_DESCRIBE_MESH_APP
#define UTOPIA_DESCRIBE_MESH_APP

#include <memory>
#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {

    class DescribeMeshApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "-describe_mesh"; }

        DescribeMeshApp();
        ~DescribeMeshApp();
    };
}  // namespace utopia

#endif  // UTOPIA_DESCRIBE_MESH_APP
