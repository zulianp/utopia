#ifndef UTOPIA_DIFFERENCE_APP
#define UTOPIA_DIFFERENCE_APP

#include "utopia_FEApp.hpp"
#include <string>
#include <memory>

namespace utopia {

    class DifferenceApp final : public FEApp {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "-diff";
        }

        DifferenceApp();
        ~DifferenceApp();

    private:
        class Impl;
    };
}


#endif //UTOPIA_DIFFERENCE_APP
