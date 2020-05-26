#ifndef UTOPIA_DIFFERENCE_APP
#define UTOPIA_DIFFERENCE_APP

#include <memory>
#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {

    class DifferenceApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command() { return "-diff"; }

        DifferenceApp();
        ~DifferenceApp();

    private:
        class Impl;
    };
}  // namespace utopia

#endif  // UTOPIA_DIFFERENCE_APP
