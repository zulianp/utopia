#ifndef UTOPIA_SMOOTH_APP
#define UTOPIA_SMOOTH_APP

#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {

    class SmoothApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "-smooth"; }
    };
}  // namespace utopia

#endif  // UTOPIA_SMOOTH_APP
